PROGRAM CREATEVERTICALRESOLUTION.EXE

Task performed:

Calculates vertical resolution for 3D modelling based on a vertical profile for a 
property, by considering the same property change with depth for all layers.

Allows the user to choose between providing the desired layer number or the 
desired property change.  


Main features:

- interaction with input file is performed through a nomfich file (input file 
indicated in IN_MODEL; reference path is indicated in ROOT_SRT);

- input file:

	- path to the file with property profile:

		PROFILEFILENAME : ...

	- path to the output geometry file:

		OUTPUTFILENAME : ...

	- vertical domain type intended (not mandatory):

		DOMAIN_TYPE : ... (only current option and default is CARTESIAN)

	- maximum bathymetry depth:

		MAXIMUMDEPTH : ...

	- required property change per layer (specified instead of LAYERS keyword)

		PROPERTY_STEP : ...
	
	- required number of layers (specified instead of PROPERTY_STEP keyword)

		LAYERS : ...

		(mandatory to provide one of LAYERS or PROPERTY_STEP)

	- maximum depth for surface layer (not mandatory, used only if PROPERTY_STEP):

		MAXDEPTH_SURFACE : ... (default is 100.0m)

	- maximum depth for bottom layer (not mandatory, used only if PROPERTY_STEP):

		MAXDEPTH_BOTTOM : ... (default is 500.0m)


- profile file:

	- block of profile depths (from lower depth to higher depth):

		<BeginDepth>
		...
		<EndDepth>

	- block of profile values (from lower depth to higher depth):

		<BeginProfileValues>
		...
		<EndProfileValues>

Some notes:

- in case the surface layer exceeds the maximum depth for surface layer for PROPERTY_STEP option
  the layer is divided till sublayers have depth smaller than maximum depth; the surface 
  layer depth is for the check added to the first profile depth (the one most near surface);

- if maximum depth of bathymetry exceeds the deeper profile depth then exceeding depth is 
  considered in the last layer (the one most in bottom); 

- in case the bottom layer exceeds the maximum depth for bottom layer for PROPERTY_STEP option
  the layer is divided till sublayers have depth smaller than maximum depth.






