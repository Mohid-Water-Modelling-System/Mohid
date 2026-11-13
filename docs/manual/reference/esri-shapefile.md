# ESRI Shapefile

## How to change the shapefile projection?

It can be done with ArcMap.

1. Add layer
2. Data Management Tools / Projections and Transforms / Feature / Project
3. Right click layer / Data / Export Data

## In CORINE shapefiles, why can't the field Code\_00 be selected as color scale in MOHID GIS?

Because it's text! It must be converted to integer, by adding a new field to the attributes table. In ArcMap:

1. Right click layer / Open Attribute Table / Add Field /Short Integer (e.g. Code\_00i)
2. Open the dbf file with Office 2003 or Access 2007 (office 2007 doesn't allow to save as dbf IV)
3. Copy the field to be saved as integer.

## How to clip with a rectangle?
