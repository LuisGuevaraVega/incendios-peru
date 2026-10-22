"""Área quemada con Sentinel-2 L2A (Microsoft Planetary Computer, sin clave).

NBR  = (NIR − SWIR2) / (NIR + SWIR2)       con B8A (865 nm) y B12 (2 190 nm), 20 m
dNBR = NBR_antes − NBR_después             ≥ 0,10 = quemado, severidad baja (Key y Benson, USGS)

Tres decisiones que cambian las hectáreas y por eso se reportan siempre:
1. Fechas antes/después y su % de nubes dentro de la zona (no el de la escena entera)
2. Máscara de nubes con la capa SCL
3. «Sembrado»: solo cuentan las manchas quemadas que tocan un foco VIIRS del incendio
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import odc.stac
import pandas as pd
import planetary_computer as pc
import pystac_client
import rasterio.features
import xarray as xr
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent))
from config import (CLASES_SEVERIDAD, CRS_METRICO, DISTANCIA_CONTROL_M,  # noqa: E402
                    MANCHA_MINIMA_HA, PERCENTIL_RUIDO, SCL_MALOS, UMBRAL_QUEMADO)

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
RES = 20               # m: la resolución de B8A, B12 y SCL
HA_PIXEL = RES * RES / 1e4
RADIO_SEMILLA_M = 750  # medio píxel VIIRS (375 m) más el error de posición


def catalogo():
    return pystac_client.Client.open(STAC, modifier=pc.sign_inplace)


def cargar(coleccion: str, bbox, fecha, bandas, res=RES, crs=f"EPSG:{CRS_METRICO}", **kw) -> xr.Dataset:
    """Carga una fecha (o rango) en la grilla común: misma bbox, CRS y resolución."""
    items = catalogo().search(collections=[coleccion], bbox=bbox, datetime=fecha).item_collection()
    if not items:
        raise ValueError(f"sin imágenes {coleccion} para {fecha}")
    ds = odc.stac.load(items, bands=bandas, bbox=bbox, crs=crs, resolution=res,
                       groupby="solar_day", chunks={}, **kw).compute()
    ds.attrs["items"] = [i.id for i in items]
    ds.attrs["baseline"] = [i.properties.get("s2:processing_baseline", "") for i in items]
    return ds


def reflectancia(da: xr.DataArray, baseline: list) -> xr.DataArray:
    """DN → reflectancia. Desde la línea base 04.00 (ene-2022) el DN trae +1000 de desfase."""
    desfase = 1000 if any(float(b or 0) >= 4 for b in baseline) else 0
    r = (da.astype("float32") - desfase) / 10000
    return r.where(da > 0)


def escena(bbox, fecha: str) -> dict:
    """NBR de una fecha con nubes enmascaradas, y el % de nubes en la zona."""
    ds = cargar("sentinel-2-l2a", bbox, fecha, ["B8A", "B12", "SCL"]).isel(time=0)
    nir, swir = (reflectancia(ds[b], ds.attrs["baseline"]) for b in ("B8A", "B12"))
    scl = ds.SCL.values
    limpio = ~np.isin(scl, SCL_MALOS) & np.isfinite(nir.values) & np.isfinite(swir.values)
    nbr = ((nir - swir) / (nir + swir)).where(limpio)
    con_dato = scl > 0
    return {"fecha": fecha, "nbr": nbr, "limpio": limpio, "agua": scl == 6,
            "nubes_pct": float(100 * (1 - limpio[con_dato].mean())) if con_dato.any() else np.nan,
            "items": ds.attrs["items"]}


def compuesto(bbox, fechas: list) -> dict:
    """Para cada píxel, la primera fecha limpia de la lista (ordenada de más a menos preferida)."""
    escenas = [escena(bbox, f) for f in fechas]
    nbr = escenas[0]["nbr"].copy()
    origen = np.where(escenas[0]["limpio"], 0, -1)
    for k, e in enumerate(escenas[1:], start=1):
        hueco = origen < 0
        nbr = nbr.where(~(hueco & e["limpio"]), e["nbr"])
        origen = np.where(hueco & e["limpio"], k, origen)
    # agua en cualquier fecha usada: los ríos cambian de nivel y simulan «quemado»
    agua = np.zeros(origen.shape, bool)
    for k, e in enumerate(escenas):
        agua |= e["agua"] & (origen == k)
    return {"nbr": nbr, "origen": origen, "escenas": escenas, "agua": agua,
            "cobertura_pct": float(100 * (origen >= 0).mean())}


def mascara_focos(plantilla: xr.DataArray, focos: gpd.GeoDataFrame, radio=RADIO_SEMILLA_M) -> np.ndarray:
    """Píxeles a ≤ radio de algún foco VIIRS (las «semillas»)."""
    buf = focos.to_crs(CRS_METRICO).buffer(radio)
    return rasterio.features.rasterize(
        [(g, 1) for g in buf], out_shape=plantilla.shape,
        transform=plantilla.odc.geobox.affine, fill=0).astype(bool)


def suavizar(dnbr: xr.DataArray) -> np.ndarray:
    """Mediana 3×3: con umbral bajo, sin esto el ruido de píxel a píxel suma hectáreas."""
    v = dnbr.values
    lleno = np.where(np.isfinite(v), v, 0.0)  # hueco = «sin cambio», lo conservador
    return np.where(np.isfinite(v), ndimage.median_filter(lleno, size=3), np.nan)


def area_quemada(dnbr: xr.DataArray, focos: gpd.GeoDataFrame, umbral=UMBRAL_QUEMADO, sembrar=True):
    """Máscara quemada: dNBR ≥ umbral, manchas ≥ 1 ha y, si se siembra, que toquen un foco."""
    quemado = suavizar(dnbr) >= umbral
    etiquetas, n = ndimage.label(quemado, structure=np.ones((3, 3)))
    tamanos = np.bincount(etiquetas.ravel(), minlength=n + 1) * HA_PIXEL
    validas = tamanos >= MANCHA_MINIMA_HA
    validas[0] = False
    if sembrar:
        semillas = mascara_focos(dnbr, focos)
        tocadas = np.zeros(n + 1, bool)
        tocadas[np.unique(etiquetas[semillas & quemado])] = True
        validas &= tocadas
    return validas[etiquetas]


def clase_severidad(dnbr: xr.DataArray, quemado: np.ndarray) -> np.ndarray:
    """0 = no quemado; 1..4 = clases de CLASES_SEVERIDAD (sobre el dNBR suavizado)."""
    v = suavizar(dnbr)
    clase = np.zeros(v.shape, "uint8")
    for i, (lo, hi, _) in enumerate(CLASES_SEVERIDAD, start=1):
        clase[(v >= lo) & (v < hi) & quemado] = i
    # un píxel quemado bajo el umbral (por la mediana) entra a la clase más baja
    clase[quemado & (clase == 0)] = 1
    return clase


def severidad(dnbr: xr.DataArray, quemado: np.ndarray) -> pd.DataFrame:
    clase = clase_severidad(dnbr, quemado)
    return pd.DataFrame([{"clase": nombre, "dnbr_desde": lo, "ha": (clase == i).sum() * HA_PIXEL}
                         for i, (lo, _, nombre) in enumerate(CLASES_SEVERIDAD, start=1)])


def perimetro(mascara: np.ndarray, plantilla: xr.DataArray, min_ha=1.0) -> gpd.GeoDataFrame:
    """Vectoriza la máscara quemada (descarta manchas menores a min_ha)."""
    formas = rasterio.features.shapes(mascara.astype("uint8"), mask=mascara,
                                      transform=plantilla.odc.geobox.affine)
    g = gpd.GeoDataFrame(geometry=[__import__("shapely").geometry.shape(s) for s, _ in formas],
                         crs=CRS_METRICO)
    g = g[g.area >= min_ha * 1e4]
    return g


def dnbr_par(pre: dict, post: dict, excluir: np.ndarray | None = None) -> xr.DataArray:
    """dNBR sin agua (SCL) ni píxeles sin combustible (p. ej. playas de río, roca)."""
    d = (pre["nbr"] - post["nbr"]).rename("dnbr")
    fuera = pre["agua"] | post["agua"]
    if excluir is not None:
        fuera = fuera | excluir
    return d.where(~fuera)


def ruido(dnbr: xr.DataArray, focos: gpd.GeoDataFrame) -> dict:
    """dNBR en zonas de control (sin focos a > 3 km): cuánto «quema» el puro ruido."""
    v = suavizar(dnbr)
    control = ~mascara_focos(dnbr, focos, DISTANCIA_CONTROL_M) & np.isfinite(v)
    c = v[control]
    return {"pixeles_control_pct": float(100 * control.mean()), "media": float(c.mean()),
            "p95": float(np.percentile(c, PERCENTIL_RUIDO)),
            "umbral": float(max(UMBRAL_QUEMADO, np.percentile(c, PERCENTIL_RUIDO)))}


def focos_entre(focos: gpd.GeoDataFrame, pre: list, post: list) -> gpd.GeoDataFrame:
    """Focos entre la pasada de Sentinel-2 (~10:30 h) de la imagen previa y la posterior.

    Un foco anterior a la imagen previa ya quemó antes: no puede sembrar este dNBR.
    """
    pasada = pd.Timedelta(hours=10, minutes=30)
    desde, hasta = pd.Timestamp(max(pre)) + pasada, pd.Timestamp(max(post)) + pasada
    return focos[focos.hora_lima.between(desde, hasta)]


def analizar(cfg: dict, focos: gpd.GeoDataFrame, excluir: np.ndarray | None = None) -> dict:
    """Corre el flujo completo de un incendio y devuelve todo lo necesario para exportar.

    excluir: máscara de píxeles que no pueden arder (en la misma grilla), o None.
    """
    bbox = cfg["bbox"]
    pre = compuesto(bbox, cfg["pre"])
    post = compuesto(bbox, cfg["post"])
    dnbr = dnbr_par(pre, post, excluir)
    r = ruido(dnbr, focos)
    semillas = focos_entre(focos, cfg["pre"], cfg["post"])
    quemado = area_quemada(dnbr, semillas, r["umbral"])
    return {"pre": pre, "post": post, "dnbr": dnbr, "quemado": quemado, "ruido": r, "excluir": excluir,
            "umbral": r["umbral"], "clase": clase_severidad(dnbr, quemado),
            "severidad": severidad(dnbr, quemado), "focos_semilla": semillas}


def sensibilidad(cfg: dict, focos: gpd.GeoDataFrame, base: dict) -> pd.DataFrame:
    """La lección de Curahuasi: el mismo incendio con otras decisiones razonables."""
    filas = []
    dnbr = base["dnbr"]
    for umbral in sorted({0.10, round(base["umbral"], 2), 0.27, 0.44}):
        for sembrar in (False, True):
            m = area_quemada(dnbr, base["focos_semilla"], umbral, sembrar)
            filas.append({"posterior": " + ".join(cfg["post"]), "umbral": umbral,
                          "sembrado_focos": sembrar, "ha": m.sum() * HA_PIXEL})
    alt = compuesto(cfg["bbox"], cfg["post_alternativa"])
    dnbr_alt = dnbr_par(base["pre"], alt, base["excluir"])
    semillas_alt = focos_entre(focos, cfg["pre"], cfg["post_alternativa"])
    for sembrar in (False, True):
        m = area_quemada(dnbr_alt, semillas_alt, base["umbral"], sembrar)
        filas.append({"posterior": " + ".join(cfg["post_alternativa"]), "umbral": round(base["umbral"], 2),
                      "sembrado_focos": sembrar, "ha": m.sum() * HA_PIXEL,
                      "nubes_pct": alt["escenas"][0]["nubes_pct"]})
    return pd.DataFrame(filas)


def rgb_compuesto(bbox, fechas: list, bandas: list, crs: str, res: float) -> np.ndarray:
    """Imagen de 3 bandas para mirar (no para medir): cada píxel de la primera fecha limpia.

    Lo que sigue nublado en todas las fechas se rellena con la primera, para no dejar huecos.
    """
    salida, hueco = None, None
    for fecha in fechas:
        ds = cargar("sentinel-2-l2a", bbox, fecha, bandas + ["SCL"], res=res, crs=crs,
                    resampling={"SCL": "nearest", "*": "bilinear"}).isel(time=0)
        a = np.dstack([reflectancia(ds[b], ds.attrs["baseline"]).values for b in bandas])
        limpio = ~np.isin(ds.SCL.values, SCL_MALOS)
        if salida is None:
            salida, hueco = a.copy(), ~limpio
            continue
        usar = hueco & limpio
        salida[usar] = a[usar]
        hueco &= ~limpio
    return salida


def nubes_por_fecha(bbox, desde: str, hasta: str) -> pd.DataFrame:
    """% de nubes dentro de la zona para cada pasada (el % de la escena entera engaña)."""
    items = catalogo().search(collections=["sentinel-2-l2a"], bbox=bbox,
                              datetime=f"{desde}/{hasta}").item_collection()
    ds = odc.stac.load(items, bands=["SCL"], bbox=bbox, crs=f"EPSG:{CRS_METRICO}", resolution=60,
                       groupby="solar_day", chunks={}).compute()
    escena = {str(i.datetime.date()): i.properties.get("eo:cloud_cover") for i in items}
    filas = []
    for t in range(ds.sizes["time"]):
        scl = ds.SCL.isel(time=t).values
        con_dato = scl > 0
        fecha = str(ds.time.values[t])[:10]
        filas.append({"fecha": fecha,
                      "nubes_zona_%": 100 * np.isin(scl, SCL_MALOS)[con_dato].mean() if con_dato.any() else np.nan,
                      "nubes_escena_%": escena.get(fecha)})
    return pd.DataFrame(filas).round(1)
