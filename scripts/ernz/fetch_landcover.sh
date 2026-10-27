#!/bin/bash
# Land Cover 2021 (LIS-L) over the Ernz bbox from the INSPIRE WFS, in pages of 10000 features, each verified and converted to GPKG.
# usage: fetch_landcover.sh <out dir> <bbox xmin,ymin,xmax,ymax in EPSG:2169> [start offset] [end offset]  (pages already converted are skipped)
O=$1; BB=$2; URL="https://wms.inspire.geoportail.lu/geoserver/lc/wfs?service=WFS&version=2.0.0&request=GetFeature&typeNames=lc:LC.LISL-LandCover2021&outputFormat=application/json&srsName=EPSG:2169&count=10000"
GDAL=/opt/homebrew/Caskroom/miniconda/base/bin
N=$(curl -s "${URL}&resultType=hits&bbox=$BB,urn:ogc:def:crs:EPSG::2169" | grep -o 'numberMatched="[0-9]*"' | grep -o '[0-9]*'); echo "features $N"
S0=${3:-0}; S1=${4:-$N}
for ((s=S0; s<S1 && s<N; s+=10000)); do
  [ -f $O/lc_$s.gpkg ] && { echo "page $s exists"; continue; }
  f=$O/lc_$s.json
  for try in 1 2 3 4; do
    curl -s --retry 2 "${URL}&startIndex=$s&bbox=$BB,urn:ogc:def:crs:EPSG::2169" -o $f
    if python3 -c "import json,sys; d=json.load(open('$f')); sys.exit(0 if len(d['features'])>0 else 1)" 2>/dev/null; then break; fi
    echo "page $s retry $try"; sleep 5
  done
  OGR_GEOJSON_MAX_OBJ_SIZE=0 $GDAL/ogr2ogr -f GPKG $O/lc_$s.gpkg $f -nln lc 2>/dev/null && rm -f $f && echo "page $s ok" || echo "page $s FAILED"
done
