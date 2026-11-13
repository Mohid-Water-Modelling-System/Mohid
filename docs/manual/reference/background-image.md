# Background Image

## What it is

A background image data item consists of an image and its geo-reference, given in the \*.img file:

```
IMAGE_FILE                    : C:\temp\389_.jpg  --> Path to the image file
X                             : 0                 --> Coordinate X of the Lower Left Corner
Y                             : 0                 --> Coordinate Y of the Lower Left Corner
WIDTH                         : 16008             --> Width in distance units
HEIGHT                        : 10011             --> Height in distance units
MIP_MAPPED                    : 1
MAINTAIN_RATIO                : 1
ALPHA_VALUE                   : 1
```

## Import ESRI geo-referenced image in MOHID GIS

### Automatic Method

ESRI geo-referenced images are comprised of two files: the image (e.g. image.jpg) and the "world file" with the coordinates (e.g. image.jpgw):

```
3                       --> pixels/unit distance X
0.000000000000
0.000000000000
-3                      --> pixels/unit distance Y
-111989.9               --> Coordinate X of Upper Left Corner
-70266.15               --> Coordinate Y of Upper Left Corner
```

This are the actions to import geo-referenced images to MOHID GIS

1. Add New Item.
2. Background Image; choose file name.
3. In this new layer created, if a window does not appear automatically, go to Properties (F4, or right-click, or Menu Data Items/Properties).
4. Choose the image file (browse)
5. Click "Import ArcView GeoRef" to import automatically the "world file"
6. Data Items -> Save All, in order to save the just created .img file

Now, everytime you want to open that background image go to Open->Background Image (type of file) and select the .img file created.

**Note:** If loading the file fails and all that you see is a white image (but in the right place) the problem is the picture resolution is too high for Mohid Gis to load it. Reduce its resolution by 50% and change the jpgw doubling pixel/unit distance X and pixel/unit distance X.

### Manual Method

Alternatively you can convert by hand the world file to Mohid GIS img file,

- WIDTH (HEIGHT) = number of width (height) pixels time pixels/unit distance X (Y). The pixel information of the image can be seen in [GIMP](http://www.gimp.org/), for example.
- The Y coordinate must be decreased by *height pixels times pixels/unit distance Y* to go from Upper Left Corner (ESRI) to Lower Left Corner (Mohid GIS)
- Write the path to the file in IMAGE\_FILE
- open the .img file created in MOHID GIS as a background image (file type).
