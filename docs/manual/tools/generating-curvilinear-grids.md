# Generating curvilinear grids

## Add a [Polygon](../io/formats/polygons.md)

Add to your [MOHID GIS](gui/mohid-gis.md) project the polygon that will be the boundary of your meshing domain

## Render Vertices

In order to 'prepare' your polygon for meshing you will need to define corners (this will be explained later). So you need to activate 'render vertices' to manipulate the vertices of the polygon. For that, right-click your polygon in the dataitem list, choose 'properties' and then check the 'render vertices' checkbox. Verify that the size of the font associated to the polygon is in agreement with the scales you are using (otherwise you might see nothing).

## Define Corners

Next step is to define the corners of the polygon. Corners are used to pave a domain in simpler sub-domains for meshing.

You always have to mark at least 4 right turn corners and this corresponds to a simple polygon. For a complex polygon one should decompose the domain in many sub-domains.The topology of any meshed sub-domain is a four edges polygon and it is defined with 2 left turns and 2 right turns (cf. picture: upper domain on the second example). Your polygon should have at least 4 right turns and the sum of right turns minus the sum of left turns should always be 4. Otherwise the program gives an error.

To define corners go to "Action->Polygons->Corners Definition". One left click on a vertex turns it into a right turn. A 'R' should appear near the vertex - if it doesn't, check size of font in polygon properties as explained before. A second click will turn it into a left turn ('R' turn to 'L'). A third click will turn the vertex to a non-corner vertex (the letter dissapears). An so on.

http://xmages.net/storage/10/1/0/9/1/upload/a772687b.png

## Generate the grid

Go to "DataItems->New". Check the "Curvilinear Grid" option and OK. Choose values for NX, NY and NPPE and then OK. NX and NY are the number of intervals of your grid in the X and Y direction respectively. NPPE is the Number of Points Per Edge. It is the minimum number of grid vertices you will have between two consecutive vertices of your boundary.

## External references

- [mohid forum topic by Hernani](http://www.mohid.com/forum/viewtopic.php?f=2&t=26)
