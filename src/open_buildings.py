"""Google Open Buildings v3: recorta las edificaciones de cada zona de incendio.

Google publica un CSV comprimido por celda S2 de nivel 4 (cientos de MB). No hace falta
guardarlo entero: se lee por partes desde la URL y solo se conservan las filas
dentro de la caja de cada incendio.

Uso:  uv run python src/open_buildings.py
"""
import sys

import geopandas as gpd
import pandas as pd
from shapely import wkt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from config import INCENDIOS, INTERIM  # noqa: E402

URL = "https://storage.googleapis.com/open-buildings-data/v3/polygons_s2_level_4_gzip/{token}_buildings.csv.gz"


def recortar(clave: str, cfg: dict) -> gpd.GeoDataFrame:
    destino = INTERIM / f"edificios_{clave}.parquet"
    if destino.exists():
        return gpd.read_parquet(destino)
    o, s, e, n = cfg["bbox_edificios"]
    partes = []
    lector = pd.read_csv(URL.format(token=cfg["s2_token"]), compression="gzip", chunksize=500_000)
    for i, bloque in enumerate(lector):
        dentro = bloque[bloque.longitude.between(o, e) & bloque.latitude.between(s, n)]
        if len(dentro):
            partes.append(dentro)
        print(f"  {clave}: bloque {i + 1} leído, {sum(map(len, partes)):,} edificaciones en la zona", flush=True)
    df = pd.concat(partes, ignore_index=True)
    gdf = gpd.GeoDataFrame(
        df.drop(columns=["geometry"]),
        geometry=df["geometry"].map(wkt.loads),
        crs=4326,
    )
    gdf.to_parquet(destino)
    return gdf


if __name__ == "__main__":
    for clave, cfg in INCENDIOS.items():
        g = recortar(clave, cfg)
        print(f"{clave}: {len(g):,} edificaciones guardadas")
