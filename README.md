# Harmony

Python requirements in requirements.txt: pip install -r requirements.txt

GTFS data for Athens came from OpenMobilityData:
https://transitfeeds.com/l/431-greece

Requires data files for Athens:
zonesDataFilename = 'athens/Data_LUTI_Athens/zones_data.csv'
zoneCentroidsShapefilenameWGS84 = 'athens/Zone_boundaries/zone_centroids_WGS84.shp'
zoneBoundariesShapefilenameWGS84 = 'athens/Zone_boundaries/zones_boundaries_WGS84.shp'

GTFS and public transport network processing files in the QUANT directory come from the QUANT project.

Run with:

python main.py
