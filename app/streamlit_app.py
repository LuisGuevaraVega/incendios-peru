"""Dashboard Streamlit de los incendios de Chichas y Río Tambo.

Uso:  uv run streamlit run app/streamlit_app.py
Lee lo que exporta src/pipeline.py (data/processed y docs/data); no recalcula nada.
"""
import json
from pathlib import Path

import folium
import geopandas as gpd
import pandas as pd
import streamlit as st
from streamlit_folium import st_folium

RAIZ = Path(__file__).resolve().parents[1]
PROC, WEB, FIGS = RAIZ / "data" / "processed", RAIZ / "docs" / "data", RAIZ / "figs"
COLORES = {"Dentro del área quemada": "#111111", "Alta": "#d7191c", "Media": "#f08a4b",
           "Baja": "#f6d55c", "Fuera de 2 km": "#9aa7b0"}

st.set_page_config(page_title="Incendios Chichas y Río Tambo", page_icon="🔥", layout="wide")
st.title("Chichas y Río Tambo, vistos desde el espacio")
st.caption("Datos al 20 de octubre de 2026 · Sentinel-2, WorldCover, Open-Meteo y Google Open Buildings · "
           "Luis Gabriel Guevara Vega")

clave = st.sidebar.radio("Incendio", ["chichas", "rio_tambo"], format_func=lambda k: {"chichas": "Chichas (Arequipa)", "rio_tambo": "Río Tambo (Junín)"}[k])
r = json.loads((PROC / clave / "resumen.json").read_text(encoding="utf-8"))
st.sidebar.markdown(f"**Imagen previa:** {', '.join(r['pre'])}  \n**Posterior:** {' → '.join(r['post'])}  \n"
                    f"**Nubes después:** {r['nubes_post_pct']} %  \n**Umbral dNBR:** {r['umbral_dnbr']}")
st.sidebar.markdown("[Dashboard 3D](https://claude.ai/artifact/5zC59rYbPi2SQ47PrArLYL)")

c = st.columns(4)
c[0].metric("Área quemada", f"{r['ha_quemadas']:,.0f} ha")
c[1].metric(r["por_cobertura"][0]["cobertura"], f"{r['por_cobertura'][0]['%']:.0f} %")
c[2].metric("Edificaciones a ≤ 1 km", f"{r['edificaciones_1km']:,}")
c[3].metric("Exposición alta", f"{r['edificaciones']['Alta']:,}")

izq, der = st.columns([3, 2])
with izq:
    st.subheader("Mapa")
    o, s, e, n = r["bbox"]
    m = folium.Map(location=[(s + n) / 2, (o + e) / 2], zoom_start=13, tiles="Esri.WorldImagery")
    folium.raster_layers.ImageOverlay(str(WEB / clave / "severidad.png"), bounds=[[s, o], [n, e]], opacity=0.75,
                                      name="Severidad").add_to(m)
    folium.GeoJson(json.loads((PROC / clave / "perimetro_quemado.geojson").read_text(encoding="utf-8")),
                   style_function=lambda _: {"color": "white", "weight": 1.5, "fillOpacity": 0}, name="Perímetro").add_to(m)
    ed = gpd.read_file(PROC / clave / "edificaciones_expuestas.geojson")
    ed = ed[ed.categoria != "Fuera de 2 km"]
    pt = ed.representative_point()
    for (_, f), p in zip(ed.iterrows(), pt):
        folium.CircleMarker([p.y, p.x], radius=4, color=COLORES[f.categoria], fill=True, fill_opacity=0.9,
                            tooltip=f"{f.categoria} · índice {f.indice_exposicion} · {f.dist_fuego_m:,.0f} m · {f.full_plus_code}").add_to(m)
    folium.LayerControl().add_to(m)
    st_folium(m, height=560, use_container_width=True)
with der:
    st.subheader("¿Qué se quemó?")
    st.dataframe(pd.DataFrame(r["por_cobertura"]), hide_index=True, use_container_width=True)
    st.subheader("Por distrito")
    st.dataframe(pd.DataFrame(r["por_distrito"]), hide_index=True, use_container_width=True)
    st.subheader("Severidad")
    st.bar_chart(pd.DataFrame(r["severidad"]).set_index("clase")["ha"])

st.subheader("Edificaciones más expuestas (Google Open Buildings)")
st.caption("Índice 0–100: cercanía al fuego (50), a favor del viento (20), cuesta arriba (15), vegetación (15). "
           "Ordena dónde mirar; no predice daños. Copia el Plus Code en Google Maps.")
tabla = pd.read_csv(PROC / clave / "edificaciones_expuestas.csv")
cols = ["full_plus_code", "categoria", "indice_exposicion", "dist_fuego_m", "a_favor_del_viento", "cuesta_arriba",
        "vegetacion_100m", "area_in_meters"]
st.dataframe(tabla[tabla.categoria != "Fuera de 2 km"].sort_values("indice_exposicion", ascending=False)[cols],
             hide_index=True, use_container_width=True)

st.subheader("La decisión importa: sensibilidad de las hectáreas")
st.dataframe(pd.DataFrame(r["sensibilidad"]), hide_index=True, use_container_width=True)
st.image(str(FIGS / "clima.png"), caption="Clima diario y pronóstico (Open-Meteo)")
