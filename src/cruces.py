"""Cruces del área quemada con otros datos abiertos.

- ESA WorldCover 2021 (10 m)       → ¿qué se quemó? bosque, pastizal, cultivo
- Copernicus DEM GLO-30            → altitud y pendiente (y el relieve 3D)
- Límites distritales INEI         → hectáreas por distrito
- Open-Meteo (sin clave)           → humedad, viento y lluvia: ¿puede crecer?
- Google Open Buildings v3         → ¿qué edificaciones están expuestas?

La exposición de una edificación NO es una predicción de que se va a quemar: es un
índice relativo (0–100) que ordena dónde mirar primero, con reglas explícitas.
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.features
import requests
import xarray as xr
from scipy import ndimage
from shapely.ops import nearest_points, unary_union

sys.path.insert(0, str(Path(__file__).parent))
from config import CORTE, CRS_METRICO, INTERIM, RAW  # noqa: E402
from sentinel import HA_PIXEL, RES, cargar, perimetro  # noqa: E402

WORLDCOVER = {10: "Bosque", 20: "Matorral", 30: "Pastizal", 40: "Cultivo", 50: "Construido",
              60: "Suelo desnudo o escaso", 70: "Nieve", 80: "Agua", 90: "Humedal",
              95: "Manglar", 100: "Musgo y liquen"}
VEGETACION = [10, 20, 30, 40, 90]
SIN_COMBUSTIBLE = [60, 70, 80]  # suelo desnudo, nieve, agua: no arden


# ------------------------------------------------------------------ raster en la grilla común
def worldcover(bbox) -> xr.DataArray:
    ds = cargar("esa-worldcover", bbox, "2021-01-01/2021-12-31", ["map"], resampling="nearest")
    return ds["map"].isel(time=0)


def dem(bbox, res=RES, crs=None) -> xr.DataArray:
    kw = {"crs": crs} if crs else {}
    ds = cargar("cop-dem-glo-30", bbox, None, ["data"], res=res, resampling="bilinear", **kw)
    return ds["data"].isel(time=0).astype("float32")


def pendiente(z: xr.DataArray) -> np.ndarray:
    gy, gx = np.gradient(z.values, RES)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def quemado_por_cobertura(quemado: np.ndarray, cobertura: xr.DataArray) -> pd.DataFrame:
    v = cobertura.values[quemado]
    t = pd.Series(v).value_counts().rename_axis("codigo").reset_index(name="pixeles")
    t["cobertura"] = t.codigo.map(WORLDCOVER)
    t["ha"] = t.pixeles * HA_PIXEL
    t["%"] = 100 * t.ha / t.ha.sum()
    return t[["cobertura", "ha", "%"]].sort_values("ha", ascending=False).reset_index(drop=True)


def quemado_por_distrito(quemado: np.ndarray, plantilla: xr.DataArray) -> pd.DataFrame:
    d = gpd.read_file(RAW / "peru_distrital_simple.geojson").to_crs(CRS_METRICO)
    caja = gpd.GeoSeries([__import__("shapely").geometry.box(*plantilla.odc.geobox.boundingbox)],
                         crs=CRS_METRICO)
    d = d[d.intersects(caja.iloc[0])].reset_index(drop=True)
    ids = rasterio.features.rasterize([(g, i + 1) for i, g in enumerate(d.geometry)],
                                      out_shape=plantilla.shape, transform=plantilla.odc.geobox.affine)
    cuenta = pd.Series(ids[quemado]).value_counts()
    filas = [{"distrito": d.NOMBDIST[i - 1].title(), "provincia": d.NOMBPROV[i - 1].title(),
              "departamento": d.NOMBDEP[i - 1].title(), "ha": n * HA_PIXEL}
             for i, n in cuenta.items() if i > 0]
    return pd.DataFrame(filas).sort_values("ha", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------ clima
def clima(lat: float, lon: float, dias_atras=30, dias_adelante=7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Diario (30 días atrás + pronóstico) y horario del pronóstico, hora de Lima."""
    destino = INTERIM / f"clima_{lat:.3f}_{lon:.3f}.json"
    if destino.exists():
        j = pd.read_json(destino, typ="series")
    else:
        r = requests.get("https://api.open-meteo.com/v1/forecast", timeout=60, params={
            "latitude": lat, "longitude": lon, "timezone": "America/Lima",
            "past_days": dias_atras, "forecast_days": dias_adelante,
            "daily": "temperature_2m_max,relative_humidity_2m_min,precipitation_sum,"
                     "wind_speed_10m_max,wind_gusts_10m_max,wind_direction_10m_dominant",
            "hourly": "wind_speed_10m,wind_direction_10m,relative_humidity_2m",
        })
        r.raise_for_status()
        destino.write_text(r.text, encoding="utf-8")
        j = pd.Series(r.json())
    diario = pd.DataFrame(j["daily"]).rename(columns={"time": "fecha"})
    diario["fecha"] = pd.to_datetime(diario.fecha)
    diario["pronostico"] = diario.fecha > CORTE
    # día propicio para que el fuego crezca: aire seco, viento y sin lluvia
    diario["dia_critico"] = ((diario.relative_humidity_2m_min < 30) & (diario.wind_gusts_10m_max >= 30)
                             & (diario.precipitation_sum < 1))
    horario = pd.DataFrame(j["hourly"]).rename(columns={"time": "hora"})
    horario["hora"] = pd.to_datetime(horario.hora)
    return diario, horario


def viento_dominante(horario: pd.DataFrame, desde: str, dias=3) -> dict:
    """Dirección media (de dónde viene) en las horas de la tarde, ponderada por velocidad.

    Las tardes (12–18 h) son cuando el aire está más seco y los incendios avanzan.
    """
    h = horario[(horario.hora > desde) & (horario.hora <= pd.Timestamp(desde) + pd.Timedelta(days=dias))]
    h = h[h.hora.dt.hour.between(12, 18)]
    rad = np.radians(h.wind_direction_10m)
    u, v = (h.wind_speed_10m * np.sin(rad)).sum(), (h.wind_speed_10m * np.cos(rad)).sum()
    direccion = (np.degrees(np.arctan2(u, v)) + 360) % 360
    return {"desde_grados": float(direccion), "hacia_grados": float((direccion + 180) % 360),
            "velocidad_media_kmh": float(h.wind_speed_10m.mean()),
            "humedad_min_pct": float(h.relative_humidity_2m.min())}


# ------------------------------------------------------------------ edificaciones expuestas
def _rumbo(a, b) -> float:
    """Rumbo (0 = norte, sentido horario) de a hacia b, en metros proyectados."""
    return (np.degrees(np.arctan2(b.x - a.x, b.y - a.y)) + 360) % 360


def _muestrear(da: xr.DataArray, x, y):
    return da.sel(x=xr.DataArray(x), y=xr.DataArray(y), method="nearest").values


def edificaciones_expuestas(edificios: gpd.GeoDataFrame, quemado: np.ndarray, plantilla: xr.DataArray,
                            focos_recientes: gpd.GeoDataFrame, cobertura: xr.DataArray,
                            z: xr.DataArray, viento: dict) -> gpd.GeoDataFrame:
    """Índice de exposición 0–100 por edificación (no quemada) y su categoría.

    proximidad   50 pts  lineal: 50 en el borde del fuego, 0 a 2 km
    viento       20 pts  la edificación queda a ±45° de hacia dónde sopla el viento de la tarde
    pendiente    15 pts  la edificación está más arriba que el punto de fuego más cercano
    combustible  15 pts  fracción de vegetación en ~100 m alrededor (WorldCover)
    """
    e = edificios.to_crs(CRS_METRICO).copy()
    centro = e.geometry.centroid
    fuego_quemado = perimetro(quemado, plantilla, min_ha=0)
    partes = list(fuego_quemado.geometry)
    if len(focos_recientes):
        # un foco VIIRS es un píxel de 375 m: se dibuja como círculo de ~190 m
        partes += list(focos_recientes.to_crs(CRS_METRICO).buffer(190))
    fuego = unary_union(partes)
    area_quemada = unary_union(list(fuego_quemado.geometry))

    e["dentro_quemado"] = centro.within(area_quemada)
    e["dist_fuego_m"] = centro.distance(fuego)
    e["dist_quemado_m"] = centro.distance(area_quemada)
    e["dist_focos_72h_m"] = (centro.distance(unary_union(list(focos_recientes.to_crs(CRS_METRICO).geometry)))
                             if len(focos_recientes) else np.nan)

    punto_fuego = [nearest_points(fuego, c)[0] for c in centro]
    rumbo = np.array([_rumbo(p, c) for p, c in zip(punto_fuego, centro)])
    dif = np.abs((rumbo - viento["hacia_grados"] + 180) % 360 - 180)
    e["a_favor_del_viento"] = (dif <= 45) & (e.dist_fuego_m > 0)

    e["altitud_m"] = _muestrear(z, centro.x.values, centro.y.values)
    z_fuego = _muestrear(z, np.array([p.x for p in punto_fuego]), np.array([p.y for p in punto_fuego]))
    e["cuesta_arriba"] = e.altitud_m > z_fuego + 5

    veg = np.isin(cobertura.values, VEGETACION).astype("float32")
    frac = xr.DataArray(ndimage.uniform_filter(veg, size=11), coords=cobertura.coords, dims=cobertura.dims)
    e["vegetacion_100m"] = _muestrear(frac, centro.x.values, centro.y.values)

    e["indice_exposicion"] = (50 * np.clip(1 - e.dist_fuego_m / 2000, 0, 1)
                              + 20 * e.a_favor_del_viento + 15 * e.cuesta_arriba
                              + 15 * e.vegetacion_100m).round(1)
    e.loc[e.dist_fuego_m > 2000, "indice_exposicion"] = 0
    e["categoria"] = pd.cut(e.indice_exposicion, [-1, 0, 35, 60, 101],
                            labels=["Fuera de 2 km", "Baja", "Media", "Alta"]).astype(str)
    e.loc[e.dentro_quemado, "categoria"] = "Dentro del área quemada"
    e["anillo"] = pd.cut(e.dist_fuego_m, [-1, 0, 250, 500, 1000, 2000, np.inf],
                         labels=["En el fuego", "0–250 m", "250–500 m", "500 m–1 km", "1–2 km", "> 2 km"]).astype(str)
    return e.to_crs(4326)
