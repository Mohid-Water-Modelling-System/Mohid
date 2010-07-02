PROGRAM INVERTEDBAROMETER.EXE

Task performed:

Calculate mean level time series for specific locations (generally tide gauge locations)
according with the Inverted (inverse) Barometer effect: a 1 mb decrease of atmospheric 
pressure below a reference pressure results in a 1.01 cm increase in water level (e.g. 
Kantha, L., K. Whitmer, G. Born, 1994, "The Inverted Barometer Effect in Altimetry: A 
Study in the North Pacific". 
Uses HDF5 files to get time and space variable atmospheric pressure data.

Allows the transformation of gauge files to use for tidal forcing in MOHID Water,
to take into account the time evolution of the reference level of each gauge according
to the new time series produced. The gauge files thus produced are ready for MOHID
Water use.  


Main features:

- interaction with input file is performed through a nomfich file (input file 
indicated in IN_MODEL; reference path is indicated in ROOT_SRT);

- input file:

	- block of HDF5 containing surface atmospheric pressure data (may be several, 
		one for each file):

		<BeginHDF5File>
		NAME : ... (path to HDF5 file)
		<EndHDF5File> 

	- definition of time window for time serie:

		START_TIME : ...
		END_TIME : ... (yyyy mm dd hh mm ss)

	- block of gauge file (may be several, one for each file):

		<BeginGaugeFile>

		NAME : ... (-> path to gauge file)

		TIMESERIE_NAME : ... (-> path to time series without extension)
				     (the name of time serie provided is used 
				      as basis for all time series names; the name
				      of each time serie is completed in 
				      InvertedBarometer with the reference «GaugeName#»,
				      where GaugeName is the one provided in the gauge
				      block in gauge file and # is the order of the gauge
				      block in gauge file)

		WRITE_MOHID_INPUT : 0/1 (-> 1 is to create a new gauge file for MOHID
					 Water use, 0 is not to create)
				        (not mandatory, by default is 0) 

		OUTPUTNAME : ... (-> path to the gauge file created for MOHID Water 
			          use with extension; mandatory only if WRITE_MOHID_INPUT is
                                  valued 1)

		<EndGaugeFile>

	- path of grid file with HDF5 files's grid:

		GRID_FILENAME : ...

	- maximum buffer size for temporary time serie storage for each gauge location 
          (not mandatory):

		 MAX_BUFFER_SIZE : ... (an integer value, by default is 100000)

Some notes:

- reference pressure is assumed 101330.0 Pa, following Dorandeu, J., P. Le Traon, 1999,
  "Effects of Global Mean Atmospheric Pressure Variations on Mean Sea Level 
  Changes from TOPEX/Poseidon", Journal of Atmospheric and Oceanic Technology, 
  Vol. 16, No. 9, pp. 1279-1283;

- in case there are gauge locations outside the provided HDF5 grid those locations' 
  coordinates and gauge file are saved in Noinformation.dat file in the directory provided
  in ROOT_SRT in nomfich.dat; if option of writing MOHID Water input file is chosen, for 
  these gauges the original gauge block is maintained with no additional keywords;

- coordinate system in gauge files is assumed equal to coordinate system of pressure data in
  HDF5 files;

- all HDF5 data files must have the same grid, the one provided by the user;

- pressure data in HDF5 files must have rank 2D and Pa units; it shoud be in HDF5 in path
  /Results//atmospheric pressure;

- HDF5 files must have non-overlapping time periods;

- time interval between pressure data values in HDF5 files may be variable: times series are
  constructed maintaing these intervals and the time intervals between files;

- since no fix DT is assumed for pressure data, no alerts are issued for no data intervals
  between HDF5 data files;

- besides water level values each time serie have a column for covered/uncovered information,
  as required by MOHID Water; every time data is assumed covered (valued 1.0). 






