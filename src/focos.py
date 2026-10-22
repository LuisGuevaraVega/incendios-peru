"""Focos de calor VIIRS del Perú: descarga, limpieza, distritos e incendios.

Fuentes (sin clave ni pago):
- INPE Queimadas, CSV diarios de Sudamérica (guarda ~30 días)
- NASA FIRMS, CSV VIIRS «South America 7d» de los tres satélites (validación)
- Límites distritales INEI, repositorio juaneladio/peru-geojson

Uso:  uv run python src/focos.py
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from config import (CORTE, CRS_METRICO, INICIO, PROCESSED, RADIO_INCENDIO_M,  # noqa: E402
                    RAW, SATELITES_VIIRS, VENTANA_INCENDIO_H)

URL_INPE = "https://dataserver-coids.inpe.br/queimadas/queimadas/focos/csv/diario/America_Sul/focos_diario_{d}.csv"
URL_FIRMS = "https://firms.modaps.eosdis.nasa.gov/data/active_fire/{carpeta}/csv/{archivo}_South_America_7d.csv"
FIRMS = {"NPP-375": ("suomi-npp-viirs-c2", "SUOMI_VIIRS_C2"),
         "NOAA-20": ("noaa-20-viirs-c2", "J1_VIIRS_C2"),
         "NOAA-21": ("noaa-21-viirs-c2", "J2_VIIRS_C2")}
URL_GEOJSON = "https://raw.githubusercontent.com/juaneladio/peru-geojson/master/{}.geojson"


def _bajar(url: str, destino: Path, forzar: bool = False) -> Path | None:
    if destino.exists() and destino.stat().st_size > 1000 and not forzar:
        return destino
    r = requests.get(url, timeout=300)
    if r.status_code != 200:
        return None
    destino.write_bytes(r.content)
    return destino


def descargar(forzar_ultimos: int = 2) -> None:
    """Baja los CSV de INPE (UTC) del periodo y los tres CSV de FIRMS.

    Los últimos días se vuelven a bajar porque INPE los completa con pasadas tardías.
    """
    (RAW / "inpe").mkdir(exist_ok=True)
    dias = pd.date_range(INICIO, pd.Timestamp(CORTE) + pd.Timedelta(days=1))
    tareas = [(URL_INPE.format(d=dia.strftime("%Y%m%d")),
               RAW / "inpe" / f"focos_diario_{dia:%Y%m%d}.csv",
               i >= len(dias) - forzar_ultimos) for i, dia in enumerate(dias)]
    with ThreadPoolExecutor(8) as pool:
        list(pool.map(lambda t: _bajar(*t), tareas))
    for carpeta, archivo in FIRMS.values():
        _bajar(URL_FIRMS.format(carpeta=carpeta, archivo=archivo),
               RAW / f"{archivo}_South_America_7d.csv", forzar=True)
    for capa in ("peru_distrital_simple", "peru_departamental_simple"):
        _bajar(URL_GEOJSON.format(capa), RAW / f"{capa}.geojson")


def distritos() -> gpd.GeoDataFrame:
    d = gpd.read_file(RAW / "peru_distrital_simple.geojson")
    return d[["IDDIST", "NOMBDIST", "NOMBPROV", "NOMBDEP", "geometry"]]


def _a_lima(utc: pd.Series) -> pd.Series:
    return utc.dt.tz_localize("UTC").dt.tz_convert("America/Lima").dt.tz_localize(None)


def focos_inpe() -> gpd.GeoDataFrame:
    """Focos VIIRS dentro del Perú, del INICIO al CORTE en hora de Lima."""
    archivos = sorted((RAW / "inpe").glob("focos_diario_*.csv"))
    df = pd.concat([pd.read_csv(f) for f in archivos], ignore_index=True)
    df = df[df.satelite.isin(SATELITES_VIIRS)].copy()
    # un mismo foco puede venir en dos archivos diarios: se deduplica por id
    df = df.drop_duplicates("id")
    df["hora_lima"] = _a_lima(pd.to_datetime(df.data_hora_gmt))
    df["fecha"] = df.hora_lima.dt.normalize()
    df = df[df.fecha.between(INICIO, CORTE)]
    g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326)
    # el país lo decide el polígono, no la columna "pais": así no se cuelan focos de frontera
    g = gpd.sjoin(g, distritos(), how="inner", predicate="within").drop(columns="index_right")
    g["satelite"] = g.satelite.replace({"NPP-375D": "NPP-375"})
    columnas = ["id", "lat", "lon", "hora_lima", "fecha", "satelite", "frp", "numero_dias_sem_chuva",
                "risco_fogo", "IDDIST", "NOMBDIST", "NOMBPROV", "NOMBDEP", "geometry"]
    return g[columnas].reset_index(drop=True)


def focos_firms(dentro_de: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Los mismos tres satélites desde FIRMS (últimos 7 días), recortados al Perú."""
    partes = []
    for sat, (_, archivo) in FIRMS.items():
        f = pd.read_csv(RAW / f"{archivo}_South_America_7d.csv")
        f["satelite"] = sat
        partes.append(f)
    f = pd.concat(partes, ignore_index=True)
    utc = pd.to_datetime(f.acq_date + " " + f.acq_time.astype(str).str.zfill(4), format="%Y-%m-%d %H%M")
    f["hora_lima"] = _a_lima(utc)
    f["fecha"] = f.hora_lima.dt.normalize()
    f = f[f.fecha <= CORTE]
    g = gpd.GeoDataFrame(f, geometry=gpd.points_from_xy(f.longitude, f.latitude), crs=4326)
    return gpd.sjoin(g, dentro_de[["geometry"]], how="inner", predicate="within").drop(columns="index_right")


def validar(inpe: gpd.GeoDataFrame, firms: gpd.GeoDataFrame) -> pd.DataFrame:
    """Focos por día y satélite en ambas fuentes, solo en los días completos que comparten."""
    dias = sorted(set(firms.fecha.unique()) & set(inpe.fecha.unique()))[1:]  # el 1.º de FIRMS llega incompleto
    a = inpe[inpe.fecha.isin(dias)].groupby("satelite").size().rename("INPE")
    b = firms[firms.fecha.isin(dias)].groupby("satelite").size().rename("FIRMS")
    t = pd.concat([a, b], axis=1)
    t["diferencia_%"] = (t.INPE / t.FIRMS - 1) * 100
    t.attrs["dias"] = (str(dias[0])[:10], str(dias[-1])[:10])
    return t


def agrupar_incendios(focos: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Un incendio = focos encadenados a ≤ 1 km y ≤ 48 h (componentes conexas)."""
    m = focos.to_crs(CRS_METRICO)
    xy = np.c_[m.geometry.x, m.geometry.y]
    pares = cKDTree(xy).query_pairs(RADIO_INCENDIO_M, output_type="ndarray")
    t = focos.hora_lima.values
    dt_h = np.abs((t[pares[:, 0]] - t[pares[:, 1]]) / np.timedelta64(1, "h"))
    pares = pares[dt_h <= VENTANA_INCENDIO_H]
    n = len(focos)
    grafo = coo_matrix((np.ones(len(pares)), (pares[:, 0], pares[:, 1])), shape=(n, n))
    _, etiqueta = connected_components(grafo, directed=False)
    focos = focos.copy()
    focos["incendio_id"] = etiqueta
    return focos


def resumen_incendios(focos: gpd.GeoDataFrame) -> pd.DataFrame:
    r = (focos.groupby("incendio_id")
         .agg(focos=("id", "size"), inicio=("hora_lima", "min"), ultimo=("hora_lima", "max"),
              frp_total=("frp", "sum"), lat=("lat", "mean"), lon=("lon", "mean"),
              distrito=("NOMBDIST", lambda s: s.mode().iat[0]),
              provincia=("NOMBPROV", lambda s: s.mode().iat[0]),
              departamento=("NOMBDEP", lambda s: s.mode().iat[0]))
         .sort_values("focos", ascending=False))
    r["duracion_dias"] = (r.ultimo - r.inicio) / pd.Timedelta(days=1)
    return r


def main():
    descargar()
    inpe = focos_inpe()
    peru = gpd.read_file(RAW / "peru_departamental_simple.geojson")[["geometry"]]
    firms = focos_firms(peru)
    print(f"INPE VIIRS Perú: {len(inpe):,} focos del {inpe.fecha.min():%d %b} al {inpe.fecha.max():%d %b}")
    v = validar(inpe, firms)
    print(f"Validación INPE vs FIRMS ({v.attrs['dias'][0]} a {v.attrs['dias'][1]}):\n{v.round(1)}")
    focos = agrupar_incendios(inpe)
    focos.to_parquet(PROCESSED / "focos_peru.parquet")
    r = resumen_incendios(focos)
    r.to_csv(PROCESSED / "incendios_peru.csv")
    v.to_csv(PROCESSED / "validacion_inpe_firms.csv")
    print(r.head(10))


if __name__ == "__main__":
    main()
