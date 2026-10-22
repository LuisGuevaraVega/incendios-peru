"""Corre todo el análisis de los dos incendios y exporta resultados, figuras y datos web.

Uso:  uv run python src/pipeline.py            (≈ 3–5 min en CPU, la 1.ª vez más)
"""
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from odc.geo.geobox import GeoBox
from odc.geo.xr import xr_reproject
from PIL import Image
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).parent))
import cruces  # noqa: E402
import sentinel as s2  # noqa: E402
from config import (CLASES_SEVERIDAD, CORTE, CRS_METRICO, DOCS, INCENDIOS,  # noqa: E402
                    INTERIM, PROCESSED, RAW)

WEB = DOCS / "data"
RES_TEXTURA = 0.0001   # ° ≈ 11 m: textura del 3D
RES_DEM_3D = 0.0003    # ° ≈ 33 m: malla del 3D (la de Copernicus es 30 m)
COLOR_SEVERIDAD = {1: (254, 217, 118), 2: (253, 141, 60), 3: (227, 26, 28), 4: (128, 0, 38)}
CODIGO_CATEGORIA = {"Dentro del área quemada": 4, "Alta": 3, "Media": 2, "Baja": 1, "Fuera de 2 km": 0}


def procesar(clave: str, focos_peru: gpd.GeoDataFrame) -> dict:
    cfg = INCENDIOS[clave]
    o, s, e, n = cfg["bbox"]
    print(f"\n=== {cfg['nombre']} ===")
    focos = focos_peru.cx[o:e, s:n].copy()
    cobertura = cruces.worldcover(cfg["bbox"])
    ev = s2.analizar(cfg, focos, excluir=np.isin(cobertura.values, cruces.SIN_COMBUSTIBLE))
    plantilla, quemado = ev["dnbr"], ev["quemado"]
    z = cruces.dem(cfg["bbox"])
    assert cobertura.shape == plantilla.shape == z.shape, "las grillas no coinciden"

    del_incendio = ev["focos_semilla"]
    lat, lon = float(del_incendio.lat.mean()), float(del_incendio.lon.mean())
    diario, horario = cruces.clima(lat, lon)
    viento = cruces.viento_dominante(horario, f"{CORTE} 23:59")
    recientes = focos[focos.fecha >= pd.Timestamp(CORTE) - pd.Timedelta(days=2)]
    edificios = gpd.read_parquet(INTERIM / f"edificios_{clave}.parquet")
    expuestas = cruces.edificaciones_expuestas(edificios, quemado, plantilla, recientes, cobertura, z, viento)

    r = {
        "cfg": cfg, "focos": focos, "focos_incendio": del_incendio, "ev": ev, "cobertura": cobertura,
        "z": z, "diario": diario, "horario": horario, "viento": viento, "recientes": recientes,
        "expuestas": expuestas,
        "por_cobertura": cruces.quemado_por_cobertura(quemado, cobertura),
        "por_distrito": cruces.quemado_por_distrito(quemado, plantilla),
        "sensibilidad": s2.sensibilidad(cfg, focos, ev),
        "perimetro": s2.perimetro(quemado, plantilla).to_crs(4326),
    }
    guardar(clave, r)
    return r


def resumen(clave: str, r: dict) -> dict:
    ev, cfg, ex = r["ev"], r["cfg"], r["expuestas"]
    pend = cruces.pendiente(r["z"])[ev["quemado"]]
    fi = r["focos_incendio"]
    cat = ex.categoria.value_counts()
    return {
        "clave": clave, "nombre": cfg["nombre"], "distrito": cfg["distrito"].title(),
        "provincia": cfg["provincia"], "departamento": cfg["departamento"], "region": cfg["region"],
        "bbox": cfg["bbox"], "corte": CORTE,
        "pre": cfg["pre"], "post": cfg["post"],
        "nubes_pre_pct": [round(x["nubes_pct"], 1) for x in ev["pre"]["escenas"]],
        "nubes_post_pct": [round(x["nubes_pct"], 1) for x in ev["post"]["escenas"]],
        "cobertura_post_pct": round(ev["post"]["cobertura_pct"], 1),
        "umbral_dnbr": round(ev["umbral"], 3), "ruido_p95": round(ev["ruido"]["p95"], 3),
        "ha_quemadas": round(float(ev["quemado"].sum() * s2.HA_PIXEL), 1),
        "severidad": ev["severidad"].round({"ha": 1}).to_dict("records"),
        "por_cobertura": r["por_cobertura"].round(1).to_dict("records"),
        "por_distrito": r["por_distrito"].round(1).to_dict("records"),
        "sensibilidad": r["sensibilidad"].round(2).fillna("").to_dict("records"),
        "focos_incendio": int(len(fi)),
        "inicio": str(fi.hora_lima.min()), "ultimo_foco": str(fi.hora_lima.max()),
        "focos_72h": int(len(r["recientes"])),
        "frp_total_mw": round(float(fi.frp.sum()), 1),
        "altitud_quemada_m": [round(float(np.nanpercentile(r["z"].values[ev["quemado"]], p))) for p in (5, 50, 95)],
        "pendiente_mediana": round(float(np.nanmedian(pend)), 1),
        "viento": {k: round(v, 1) for k, v in r["viento"].items()},
        "dias_criticos_pronostico": int(r["diario"].query("pronostico").dia_critico.sum()),
        "edificaciones_total": int(len(ex)),
        "edificaciones": {k: int(cat.get(k, 0)) for k in CODIGO_CATEGORIA},
        "edificaciones_2km": int((ex.dist_fuego_m <= 2000).sum()),
        "edificaciones_1km": int((ex.dist_fuego_m <= 1000).sum()),
        "edificaciones_500m": int((ex.dist_fuego_m <= 500).sum()),
        "edificaciones_viento_2km": int((ex.a_favor_del_viento & (ex.dist_fuego_m <= 2000)).sum()),
    }


def guardar(clave: str, r: dict) -> None:
    d = PROCESSED / clave
    d.mkdir(exist_ok=True)
    ev = r["ev"]
    ev["dnbr"].odc.write_cog(d / "dnbr.tif", overwrite=True)
    sev = ev["dnbr"].copy(data=ev["clase"]).astype("uint8")
    sev.odc.write_cog(d / "severidad.tif", overwrite=True)
    r["perimetro"].to_file(d / "perimetro_quemado.geojson", driver="GeoJSON")
    r["focos"].to_file(d / "focos.geojson", driver="GeoJSON")
    ex = r["expuestas"]
    ex.to_file(d / "edificaciones_expuestas.geojson", driver="GeoJSON")
    ex.drop(columns="geometry").to_csv(d / "edificaciones_expuestas.csv", index=False)
    for nombre in ("por_cobertura", "por_distrito", "sensibilidad", "diario"):
        r[nombre].to_csv(d / f"{nombre}.csv", index=False)
    ev["severidad"].to_csv(d / "severidad.csv", index=False)
    (d / "resumen.json").write_text(json.dumps(resumen(clave, r), ensure_ascii=False, indent=1, default=str),
                                    encoding="utf-8")


# ------------------------------------------------------------------ datos para la web (3D)
def _png_terrarium(z: np.ndarray, destino: Path) -> None:
    """Altitud → PNG «terrarium» (R·256 + G + B/256 − 32768), el formato que lee el visor."""
    v = np.nan_to_num(z, nan=float(np.nanmin(z))) + 32768
    rgb = np.dstack([np.floor(v / 256), np.floor(v) % 256, np.floor((v - np.floor(v)) * 256)])
    Image.fromarray(rgb.astype("uint8")).save(destino)


def _sombreado(z: np.ndarray, paso_m: float, azimut=315, altura=40) -> np.ndarray:
    """Relieve sombreado (luz del noroeste) en escala de grises."""
    z = np.nan_to_num(z, nan=float(np.nanmin(z)))
    gy, gx = np.gradient(z, paso_m)
    pend = np.arctan(np.hypot(gx, gy))
    asp = np.arctan2(-gx, gy)
    az, al = np.radians(azimut), np.radians(altura)
    s = np.sin(al) * np.cos(pend) + np.cos(al) * np.sin(pend) * np.cos(az - asp)
    g = (255 * np.clip(0.15 + 0.85 * s, 0, 1)).astype("uint8")
    return np.dstack([g, g, (g * 0.97).astype("uint8")])


def _top_expuestas(ex: gpd.GeoDataFrame, n=25) -> list:
    t = ex[ex.categoria != "Fuera de 2 km"].sort_values("indice_exposicion", ascending=False).head(n)
    p = t.geometry.representative_point()
    return [{"plus_code": r.full_plus_code, "lon": round(q.x, 5), "lat": round(q.y, 5),
             "area_m2": round(r.area_in_meters), "dist_m": round(r.dist_fuego_m),
             "viento": bool(r.a_favor_del_viento), "cuesta_arriba": bool(r.cuesta_arriba),
             "vegetacion": round(float(r.vegetacion_100m), 2), "indice": float(r.indice_exposicion),
             "categoria": r.categoria} for r, q in zip(t.itertuples(), p)]


def _estirar(a: np.ndarray, lo: np.ndarray, hi: np.ndarray, gamma=0.85) -> np.ndarray:
    return (255 * np.clip((a - lo) / (hi - lo), 0, 1) ** gamma).astype("uint8")


def exportar_web(clave: str, r: dict) -> dict:
    cfg, ev = r["cfg"], r["ev"]
    bbox = cfg["bbox"]
    dest = WEB / clave
    dest.mkdir(parents=True, exist_ok=True)
    crs = "EPSG:4326"

    z = cruces.dem(bbox, res=RES_DEM_3D, crs=crs)
    _png_terrarium(z.values, dest / "dem.png")
    Image.fromarray(_sombreado(z.values, RES_DEM_3D * 111_000)).save(dest / "relieve_post.jpg", quality=88)
    caja = [float(z.longitude.min()) - RES_DEM_3D / 2, float(z.latitude.min()) - RES_DEM_3D / 2,
            float(z.longitude.max()) + RES_DEM_3D / 2, float(z.latitude.max()) + RES_DEM_3D / 2]

    texturas = {"natural": ["B04", "B03", "B02"], "infrarrojo": ["B12", "B8A", "B04"]}  # + relieve
    for nombre, bandas in texturas.items():
        post = s2.rgb_compuesto(bbox, cfg["post"], bandas, crs, RES_TEXTURA)
        pre = s2.rgb_compuesto(bbox, cfg["pre"], bandas, crs, RES_TEXTURA)
        # un solo estiramiento para las 3 bandas y para antes/después: colores comparables
        lo, hi = np.nanpercentile(pre, 2), np.nanpercentile(pre, 98)
        for cuando, a in (("post", post), ("pre", pre)):
            Image.fromarray(_estirar(np.nan_to_num(a), lo, hi)).save(dest / f"{nombre}_{cuando}.jpg", quality=86)

    # severidad reproyectada a la grilla de la textura, como capa RGBA
    sev = ev["dnbr"].copy(data=ev["clase"].astype("uint8"))
    alto, ancho = post.shape[:2]
    gb = GeoBox.from_bbox(caja, crs=crs, shape=(alto, ancho))
    sev4326 = xr_reproject(sev, gb, resampling="nearest", dst_nodata=0).values
    rgba = np.zeros((*sev4326.shape, 4), "uint8")
    for k, color in COLOR_SEVERIDAD.items():
        rgba[sev4326 == k] = (*color, 205)
    Image.fromarray(rgba).save(dest / "severidad.png")

    # focos del incendio con altitud y día; zona VIIRS acumulada por día
    f = r["focos_incendio"].copy()
    dias = pd.date_range(max(cfg["pre"]), CORTE)  # desde el día de la imagen previa
    f["dia"] = ((f.fecha - dias[0]) / pd.Timedelta(days=1)).astype(int).clip(0)
    f["alt"] = z.sel(longitude=f.lon.to_xarray(), latitude=f.lat.to_xarray(), method="nearest").values
    f["horas"] = (f.hora_lima - dias[0]) / pd.Timedelta(hours=1)
    focos = f[["lon", "lat", "alt", "horas", "frp", "dia"]].round(
        {"lon": 5, "lat": 5, "alt": 0, "horas": 1, "frp": 1}).values.tolist()
    fm = f.to_crs(CRS_METRICO)
    zonas = []
    for i in range(len(dias)):
        hasta = fm[fm.dia <= i]
        if len(hasta):
            g = gpd.GeoSeries([unary_union(list(hasta.buffer(190)))], crs=CRS_METRICO).to_crs(4326)
            g = shapely.set_precision(g.simplify(0.0006).values, 1e-5)
            zonas.append(json.loads(gpd.GeoSeries(g, crs=4326).to_json())["features"][0]["geometry"])
        else:
            zonas.append(None)
    por_dia = f.groupby("dia").agg(focos=("lon", "size"), frp=("frp", "sum")).reindex(range(len(dias)), fill_value=0)

    per = r["perimetro"].to_crs(CRS_METRICO).buffer(0).to_crs(4326).simplify(0.00015)
    per = shapely.set_precision(unary_union(list(per)), 1e-5)
    perimetro = json.loads(gpd.GeoSeries([per], crs=4326).to_json())["features"][0]["geometry"]

    ex = r["expuestas"]
    c = ex.geometry.representative_point()
    alt = z.sel(longitude=c.x.to_xarray(), latitude=c.y.to_xarray(), method="nearest").values
    edificios = [[round(x, 5), round(y, 5), float(round(a)), float(i), CODIGO_CATEGORIA[k], float(round(d))]
                 for x, y, a, i, k, d in zip(c.x, c.y, alt, ex.indice_exposicion, ex.categoria, ex.dist_fuego_m)]

    detalle = {
        **resumen(clave, r), "bounds": caja, "dem_rango": [float(np.nanmin(z)), float(np.nanmax(z))],
        "texturas": [*texturas, "relieve"], "dias": [{"fecha": str(d.date()), "focos": int(por_dia.focos[i]),
                                              "frp": round(float(por_dia.frp[i]), 1)} for i, d in enumerate(dias)],
        "focos": focos, "zonas": zonas, "perimetro": perimetro, "edificios": edificios,
        "top_expuestas": _top_expuestas(ex),
        "clima": r["diario"].assign(fecha=r["diario"].fecha.dt.strftime("%Y-%m-%d")).round(1).to_dict("records"),
        "clases_severidad": [{"clase": n, "desde": lo, "color": COLOR_SEVERIDAD[i + 1]}
                             for i, (lo, _, n) in enumerate(CLASES_SEVERIDAD)],
    }
    (dest / "detalle.json").write_text(json.dumps(detalle, ensure_ascii=False, default=str), encoding="utf-8")
    return detalle


def exportar_nacional(focos: gpd.GeoDataFrame) -> None:
    """Focos del periodo (para el mapa nacional) y serie diaria por departamento."""
    WEB.mkdir(parents=True, exist_ok=True)
    dia0 = focos.fecha.min()
    pts = np.c_[focos.lon.round(3), focos.lat.round(3), ((focos.fecha - dia0) / pd.Timedelta(days=1)).astype(int),
                focos.frp.round(0)]
    serie = focos.groupby(["fecha", "NOMBDEP"]).size().unstack(fill_value=0)
    inc = pd.read_csv(PROCESSED / "incendios_peru.csv")
    nacional = {
        "desde": str(dia0.date()), "hasta": CORTE, "total": int(len(focos)), "puntos": pts.tolist(),
        "dias": [str(d.date()) for d in serie.index], "por_dia": serie.sum(axis=1).tolist(),
        "departamentos": serie.sum().sort_values(ascending=False).to_dict(),
        "top_incendios": inc.head(10).round(2).to_dict("records"),
    }
    (WEB / "nacional.json").write_text(json.dumps(nacional, ensure_ascii=False, default=str), encoding="utf-8")
    dep = gpd.read_file(RAW / "peru_departamental_simple.geojson")[["NOMBDEP", "geometry"]]
    dep["geometry"] = dep.simplify(0.01)
    dep.to_file(WEB / "departamentos.geojson", driver="GeoJSON")


def main():
    focos = gpd.read_parquet(PROCESSED / "focos_peru.parquet")
    exportar_nacional(focos)
    for clave in INCENDIOS:
        r = procesar(clave, focos)
        d = exportar_web(clave, r)
        print(json.dumps({k: d[k] for k in ("ha_quemadas", "umbral_dnbr", "cobertura_post_pct",
                                             "edificaciones", "edificaciones_1km", "viento")},
                         ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
