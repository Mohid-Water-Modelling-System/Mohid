# Module Drawing

## Overview

This module reads several geographical type data formats that can be used in MOHID (e.g. [XYZ Points](../io/formats/xyz-points.md), [Polygons](../io/formats/polygons.md), [Lines](../io/formats/lines.md))

## Point inside a polygon algorithm

**Acceleration test**

The first action performed is an acceleration test. If the point is not inside the maximum limits of the polygon then it is not considered to be inside the polygon.

**Intersections test**

If point is inside the maximum bounds of the polygon then there is a possibility of being inside the polygon. Thus, it is drawn a semi-line on the X axis to the right of the point. If the semi-line intersects the polygon an odd number of times then the point is inside the polygon. The intersection test is performed for every segment of the polygon.

**Overlapping test**

If the point belongs to a vertix of the polygon or belongs to one of the segments which define the polygon then it is considered to be inside the polygon.

**Exceptions handling**

There are a few cases in which special care needs to be taken, regarding the intersection counting. If the intersection point coincides with one of the vertices then one intersection is counted, but only if it is the vertice with the higher Y value. This way only one intersection is counted, as one polygon vertice belongs to two segments and otherwise it would be tested twice. Another case is when the segment belongs to the semi-line and no intersection is counted.
