## DESCRIPTION

*i.up42.import* is a Python 3 script which uses the Python SDK for UP42.
UP42 is a platform for Earth data and analytics. You have to sign in to
UP42 before using this GRASS GIS addon (<https://up42.com/>). So far,
only the download of the Pléiades data is supported by this addon. The
download region can be defined by an input vector map **input** or a
GeoJSON polygon **input_geojson**. Here it is important that the area
size must be at least 0.1 sqkm. if **input_geojson** is used, the
GeoJSON must be passed as string with a polygon, UP42 only allows one
single feature and not multiple features inside a geojson file and also
no multipolygons. With the flag **p** the credits of a process are
estimated.

## EXAMPLES

### Estimate Credits for download of Pléiades data contained in a GeoJSON Polygon

```sh
# setting authentication
export UP42_PROJECT_ID=[PROJECT_ID]
export UP42_PROJECT_API_KEY=[PROJECT_API_KEY]

# polygon in GeoJSON format
GEOJSON='{"type": "FeatureCollection","features": [{"type": "Feature","properties": {},"geometry": {"type": "Polygon","coordinates": [[[7.095494270324707,50.742294306843085],[7.092854976654053,50.740447578361746],[7.098348140716552,50.73836992176511],[7.100987434387207,50.73957850446387],[7.098090648651123,50.7425794569557],[7.095494270324707,50.742294306843085]]]}}]}'

# import of Pléiades in this polygon
i.up42.import input_geojson="$GEOJSON" clouds=2 -p directory=vhr/ start=2020-04-01 end=2020-04-30 output=vhr
```

## SEE ALSO

*[r.import](https://grass.osgeo.org/grass-stable/manuals/r.import.html)*

## AUTHOR

Anika Weinmann, [mundialis](https://www.mundialis.de/), Germany
