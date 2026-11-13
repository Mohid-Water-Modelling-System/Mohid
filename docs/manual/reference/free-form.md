# Free form

A free form bathymetry uses the following keywords:

- COMMENT\_HEADER - a optional comment line describing the bathymetry;
- COMMENT\_DLX - a optional comment line;
- DLX - a array with with IMX elements, containing the length of each reservoier segment, from 1 to IMX;
- COMMENT\_ELWS - a optional comment line;
- ELWS - a array with with IMX elements, containing the inicial surface elevation at each reservoier segment, from 1 to IMX;
- COMMENT\_EPHI0 - aCOMMENT\_FRIC
- EPHI0 - a array with with IMX elements, Angle/Segment orientation, from 1 to IMX;
- COMMENT\_FRIC - a optional comment line;
- FRIC - a array with with IMX elements, containing the Bottom Friction, from 1 to IMX;
- COMMENT\_H - a optional comment line;
- H - a array with with KMX elements, Dz/Layer thickness, from 1 (surface) to KMX (bottom);
- COMMENT\_SEGMENT - a optional comment line;
- <Segments>...</Segments> - block of data with IMX raws and 4 columns. The first column specifies the segment it refers to, the second column contains the width of the surface layer, the third column contains the bottom width and the forth column specifies the bottom layer. If the last 3 columns are zero it means it is the extremity of a branch.

Sample free form bathymetry:

---

COMMENT\_HEADER : Bathymetry File for CE-Qual-W2 V3.12 - ALQUEVA (SEG 25, LAY 30, BRA 1)

COMMENT\_DLX : Segment length 

$$
DLX
$$


DLX  : 3428.67 3428.67 2142.86 2142.86 1714.29 2142.86 2142.86 1285.71 1714.29 1714.29 3000.00 1714.29 2142.86 2142.86 2785.71 3000.00 2571.43 2571.43 1714.29 2571.43 1928.57 857.14 1071.43 1928.57 1928.57

COMMENT\_ELWS : Water surface elevation 

$$
ELWS
$$


ELWS  : 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00 50.00

COMMENT\_EPHI0 : Angle/Segment orientation 

$$
PHI0
$$


PHI0  : 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00

COMMENT\_FRIC : Bottom Friction 

$$
FRICT
$$


FRIC : 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00 75.00

COMMENT\_H : Dz/Layer thickness (weigth,m) 

$$
H
$$


H : 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000 2.000

COMMENT\_SEGMENT : Dist(m) - Upstream segmento of
!i B\_K2 B\_KN KN
<Segments>
1 0 0 0
2 1500 500 3
3 1500 500 4
4 1500 500 5
5 1500 500 7
6 1500 500 8
7 1500 500 9
8 1500 500 9
9 1600 670 9
10 1600 670 10
11 1600 670 11
12 1600 670 11
13 1600 670 11
14 1600 670 11
15 1700 670 12
16 1700 670 12
17 1700 670 14
18 1700 300 14
19 1700 300 14
20 1700 300 19
21 1700 300 19
22 2000 200 25
23 2000 200 28
24 2000 200 29
25 0 0 0
</Segments>

---
