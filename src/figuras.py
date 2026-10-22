"""Figuras del README y del notebook, a partir de lo que exporta pipeline.py.

Uso:  uv run python src/figuras.py
"""
import json
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from config import CORTE, DOCS, FIGS, INCENDIOS, PROCESSED, RAW  # noqa: E402

WEB = DOCS / "data"
plt.rcParams.update({
    "figure.dpi": 110, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "font.family": "DejaVu Sans",
})
FUEGO = "#B4441A"
AZUL = "#2E5E7E"
CATEGORIAS = {  # color, tamaño, orden de dibujo
    "Fuera de 2 km": ("#9aa7b0", 4),
    "Baja": ("#F6D55C", 7),
    "Media": ("#F08A4B", 10),
    "Alta": ("#D7191C", 14),
    "Dentro del área quemada": ("#000000", 16),
}


def detalle(clave: str) -> dict:
    return json.loads((WEB / clave / "detalle.json").read_text(encoding="utf-8"))


def mapa_incendio(clave: str, ax=None, leyenda=True):
    """Imagen posterior + severidad + perímetro + edificaciones por exposición + viento."""
    d = detalle(clave)
    o, s, e, n = d["bounds"]
    base = Image.open(WEB / clave / "natural_post.jpg").convert("RGBA")
    base.alpha_composite(Image.open(WEB / clave / "severidad.png"))
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(base, extent=[o, e, s, n])
    per = gpd.read_file(PROCESSED / clave / "perimetro_quemado.geojson")
    per.boundary.plot(ax=ax, color="white", linewidth=0.7)
    ed = gpd.read_file(PROCESSED / clave / "edificaciones_expuestas.geojson")
    for cat, (color, tam) in CATEGORIAS.items():
        sub = ed[ed.categoria == cat]
        sub.representative_point().plot(ax=ax, color=color, markersize=tam, edgecolor="white", linewidth=0.3, zorder=3)
    f = gpd.read_file(PROCESSED / clave / "focos.geojson")
    f["fecha"] = pd.to_datetime(f.fecha)
    rec = f[f.fecha >= pd.Timestamp(CORTE) - pd.Timedelta(days=2)]
    rec.plot(ax=ax, marker="^", color="#00E5FF", markersize=14, edgecolor="black", linewidth=0.3, zorder=4)
    # flecha del viento de la tarde (hacia dónde sopla)
    rumbo = np.radians(d["viento"]["hacia_grados"])
    lx, ly = o + 0.12 * (e - o), n - 0.12 * (n - s)
    largo = 0.07 * (e - o)
    ax.annotate("", xy=(lx + largo * np.sin(rumbo), ly + largo * np.cos(rumbo)), xytext=(lx, ly),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=2.2))
    ax.text(lx, ly - 0.035 * (n - s), f"viento\n{d['viento']['velocidad_media_kmh']:.0f} km/h", color="white",
            ha="center", va="top", fontsize=8, fontweight="bold")
    ax.set_xlim(o, e)
    ax.set_ylim(s, n)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(False)
    ax.set_title(f"{d['nombre']} ({d['provincia']}, {d['departamento']}) · {d['ha_quemadas']:,.0f} ha quemadas",
                 fontsize=10.5, loc="left")
    if leyenda:
        h = [Patch(color=tuple(c / 255 for c in x["color"]), label=f"Severidad {x['clase'].lower()}")
             for x in d["clases_severidad"]]
        h += [Line2D([], [], marker="o", ls="", color=c, markersize=np.sqrt(t) + 1, label=f"Edificación: {k.lower()}")
              for k, (c, t) in CATEGORIAS.items()]
        h += [Line2D([], [], marker="^", ls="", color="#00E5FF", markeredgecolor="k", label="Focos 18–20 oct")]
        ax.legend(handles=h, loc="lower right", fontsize=7, framealpha=0.85)
    return ax


def antes_despues(clave: str):
    d = detalle(clave)
    o, s, e, n = d["bounds"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 5))
    for a, (archivo, titulo) in zip(ax, [
            ("natural_pre.jpg", f"Antes · {d['pre'][-1]}"),
            ("natural_post.jpg", f"Después · {' + '.join(d['post'])}"),
            ("infrarrojo_post.jpg", "Después, infrarrojo (B12-B8A-B4): quemado = rojo oscuro")]):
        a.imshow(Image.open(WEB / clave / archivo), extent=[o, e, s, n])
        a.set_title(titulo, fontsize=9, loc="left")
        a.set_xticks([])
        a.set_yticks([])
        a.grid(False)
    fig.suptitle(f"{d['nombre']}: Sentinel-2 L2A, color natural e infrarrojo", x=0.01, ha="left", fontsize=11)
    plt.tight_layout()
    return fig


def nacional():
    nac = json.loads((WEB / "nacional.json").read_text(encoding="utf-8"))
    dep = gpd.read_file(WEB / "departamentos.geojson")
    pts = np.array(nac["puntos"])
    fig = plt.figure(figsize=(12, 6.2))
    a = fig.add_axes([0.0, 0.02, 0.42, 0.92])
    dep.plot(ax=a, color="#efe8dc", edgecolor="#b8afa2", linewidth=0.4)
    a.scatter(pts[:, 0], pts[:, 1], s=0.4, c=FUEGO, alpha=0.35, linewidths=0)
    for clave, cfg in INCENDIOS.items():
        o, s, e, n = cfg["bbox"]
        a.plot([o, e, e, o, o], [s, s, n, n, s], color="black", lw=1)
        a.annotate(cfg["nombre"], ((o + e) / 2, n), xytext=(-60 if clave == "chichas" else 25, 12),
                   textcoords="offset points", fontsize=9, fontweight="bold",
                   arrowprops=dict(arrowstyle="-", color="black", lw=0.6))
    a.set_axis_off()
    a.set_title(f"{nac['total']:,} focos VIIRS en el Perú · {nac['desde']} a {nac['hasta']}", fontsize=10, loc="left")
    b = fig.add_axes([0.5, 0.56, 0.48, 0.36])
    dias = pd.to_datetime(nac["dias"])
    serie = pd.Series(nac["por_dia"], index=dias)
    b.bar(dias, serie, color=FUEGO, alpha=0.5, width=0.85)
    b.plot(dias, serie.rolling(7, min_periods=1).mean(), color=FUEGO, lw=2, label="media móvil 7 días")
    b.set_title("Focos por día (S-NPP + NOAA-20 + NOAA-21)", fontsize=9.5, loc="left")
    b.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    b.legend(frameon=False, fontsize=8)
    c = fig.add_axes([0.62, 0.07, 0.36, 0.38])
    top = pd.Series(nac["departamentos"]).head(10)[::-1]
    c.barh([x.title() for x in top.index], top.values, color=[AZUL if x in ("JUNIN", "AREQUIPA") else "#9aa7b0" for x in top.index])
    c.set_title("Departamentos con más focos", fontsize=9.5, loc="left")
    return fig


def sensibilidad():
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 6))
    for a, clave in zip(ax, INCENDIOS):
        d = detalle(clave)
        t = pd.DataFrame(d["sensibilidad"])
        etiqueta = (t.posterior.str.replace("2026-", "").str.replace(" + ", "+") + " · umbral "
                    + t.umbral.map("{:.2f}".format)
                    + "\n" + t.sembrado_focos.map({True: "sembrado con focos", False: "sin sembrar"}))
        elegido = (t.posterior == " + ".join(d["post"])) & (t.umbral.round(2) == round(d["umbral_dnbr"], 2)) & t.sembrado_focos
        a.barh(etiqueta, t.ha, color=np.where(elegido, FUEGO, "#c9c0b3"))
        for y, v in enumerate(t.ha):
            a.text(v, y, f" {v:,.0f}", va="center", fontsize=7.5)
        a.set_title(f"{d['nombre']}: el mismo incendio, {t.ha.min():,.0f} a {t.ha.max():,.0f} ha", fontsize=10, loc="left")
        a.tick_params(axis="y", labelsize=7.5)
        a.set_xlabel("hectáreas")
        a.set_xlim(0, t.ha.max() * 1.18)
    plt.tight_layout()
    return fig


def clima():
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8), sharey=False)
    for a, clave in zip(ax, INCENDIOS):
        d = detalle(clave)
        c = pd.DataFrame(d["clima"])
        c["fecha"] = pd.to_datetime(c.fecha)
        a.fill_between(c.fecha, c.relative_humidity_2m_min, color=AZUL, alpha=0.25, label="humedad mínima (%)")
        a.plot(c.fecha, c.wind_gusts_10m_max, color=FUEGO, lw=1.6, label="ráfaga máxima (km/h)")
        a.bar(c.fecha, c.precipitation_sum, color=AZUL, width=0.7, label="lluvia (mm)")
        crit = c[c.dia_critico]
        a.scatter(crit.fecha, crit.wind_gusts_10m_max, color="black", s=12, zorder=3, label="día crítico")
        a.axvspan(pd.Timestamp(CORTE) + pd.Timedelta(hours=12), c.fecha.max() + pd.Timedelta(hours=12),
                  color="#000", alpha=0.06, label="pronóstico")
        a.set_title(f"{d['nombre']}: clima diario (Open-Meteo)", fontsize=10, loc="left")
        a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        a.tick_params(axis="x", labelsize=7.5)
    ax[0].legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    plt.tight_layout()
    return fig


def exposicion():
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.6))
    orden = ["Dentro del área quemada", "Alta", "Media", "Baja"]
    for a, clave in zip(ax, INCENDIOS):
        d = detalle(clave)
        v = [d["edificaciones"][k] for k in orden]
        a.barh(orden[::-1], v[::-1], color=[CATEGORIAS[k][0] for k in orden[::-1]])
        for y, x in enumerate(v[::-1]):
            a.text(x, y, f" {x:,}", va="center", fontsize=8)
        a.set_title(f"{d['nombre']}: {d['edificaciones_2km']:,} edificaciones a ≤ 2 km del fuego "
                    f"(de {d['edificaciones_total']:,} en la zona)", fontsize=9.5, loc="left")
        a.set_xlim(0, max(v) * 1.2)
    plt.tight_layout()
    return fig


def main():
    nacional().savefig(FIGS / "nacional.png", dpi=130, bbox_inches="tight")
    sensibilidad().savefig(FIGS / "sensibilidad.png", dpi=130, bbox_inches="tight")
    clima().savefig(FIGS / "clima.png", dpi=130, bbox_inches="tight")
    exposicion().savefig(FIGS / "exposicion.png", dpi=130, bbox_inches="tight")
    for clave in INCENDIOS:
        fig, ax = plt.subplots(figsize=(8.5, 8.5))
        mapa_incendio(clave, ax)
        fig.savefig(FIGS / f"{clave}_mapa.png", dpi=130, bbox_inches="tight")
        antes_despues(clave).savefig(FIGS / f"{clave}_antes_despues.png", dpi=110, bbox_inches="tight")
    plt.close("all")
    print("figuras en", FIGS)


if __name__ == "__main__":
    main()
