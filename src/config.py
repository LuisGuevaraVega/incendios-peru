"""Configuración única del proyecto: rutas, fechas de corte y los dos incendios.

Todas las decisiones de método que cambian el resultado (fechas, umbrales, zona)
están aquí, para que se puedan reportar y cambiar en un solo lugar.
"""
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
RAW = RAIZ / "data" / "raw"
INTERIM = RAIZ / "data" / "interim"
PROCESSED = RAIZ / "data" / "processed"
FIGS = RAIZ / "figs"
DOCS = RAIZ / "docs"
for carpeta in (RAW, INTERIM, PROCESSED, FIGS):
    carpeta.mkdir(parents=True, exist_ok=True)

# Corte de los datos: hasta el 20 de octubre de 2026 inclusive (hora de Lima)
INICIO = "2026-09-21"
CORTE = "2026-10-20"

# CRS métrico para áreas y distancias (UTM 18S cubre ambos incendios)
CRS_METRICO = 32718

# Regla del profe para pasar de focos a incendios
RADIO_INCENDIO_M = 1000
VENTANA_INCENDIO_H = 48

# Satélites VIIRS (375 m). Se excluyen GOES y MODIS: otra resolución, otra física
SATELITES_VIIRS = ["NPP-375", "NPP-375D", "NOAA-20", "NOAA-21"]

# Umbral dNBR para "quemado" y clases de severidad (Key & Benson, USGS).
# Se usa 0,10 (inicio de "severidad baja") y no 0,27: en la puna y en pastizales el
# fuego deja un dNBR bajo y con 0,27 se perdía ~90 % de la cicatriz de Chichas.
# 0,27 queda en la tabla de sensibilidad. Es un piso: si el ruido del par de imágenes
# (p95 del dNBR en zonas sin focos a > 3 km) es mayor, se usa ese valor.
UMBRAL_QUEMADO = 0.10
PERCENTIL_RUIDO = 95
DISTANCIA_CONTROL_M = 3000
MANCHA_MINIMA_HA = 1.0
CLASES_SEVERIDAD = [
    (0.10, 0.27, "Baja"),
    (0.27, 0.44, "Moderada-baja"),
    (0.44, 0.66, "Moderada-alta"),
    (0.66, 9.99, "Alta"),
]
# SCL de Sentinel-2 que se enmascaran: sombra de nube, nubes media/alta, cirros, nieve
SCL_MALOS = [0, 1, 3, 8, 9, 10, 11]

# Distancias para las viviendas en riesgo (m)
ANILLOS_RIESGO = [0, 250, 500, 1000, 2000]

INCENDIOS = {
    "chichas": {
        "nombre": "Chichas",
        "distrito": "CHICHAS",
        "provincia": "Condesuyos",
        "departamento": "Arequipa",
        "region": "Sierra",
        # zona de análisis: focos del incendio + ~3 km
        "bbox": [-73.03, -15.56, -72.93, -15.43],
        "bbox_edificios": [-73.08, -15.61, -72.88, -15.38],
        "pre": ["2026-10-15"],
        # la del 20 oct es la más reciente; donde tiene nubes se completa con la del 18
        "post": ["2026-10-20", "2026-10-18"],
        # sensibilidad: solo la del 18, sin nubes pero con el fuego todavía activo
        "post_alternativa": ["2026-10-18"],
        "s2_token": "917",
    },
    "rio_tambo": {
        "nombre": "Río Tambo",
        "distrito": "RIO TAMBO",
        "provincia": "Satipo",
        "departamento": "Junín",
        "region": "Selva central",
        "bbox": [-73.94, -11.95, -73.78, -11.775],
        "bbox_edificios": [-73.99, -12.00, -73.73, -11.73],
        "pre": ["2026-09-28"],
        # selva nublada: compuesto con la última observación limpia de cada píxel
        "post": ["2026-10-20", "2026-10-18", "2026-10-13"],
        "post_alternativa": ["2026-10-13"],
        "s2_token": "911",
    },
}
