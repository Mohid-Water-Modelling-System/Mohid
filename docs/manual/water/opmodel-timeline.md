# OpModel TimeLine

This article describes the Mercator-IST operational system timeline.

## Committed Computer Resources

- Maretec-IST intranet
- Hidromod-IST intranet
- \\\\guillaume
- \\\\192.168.23.151 or \\\\Dual\_64
- ~\\\\kepler~ \\\\alexandria
- \\\\data.mohid.com

### Relevant locations

#### \\\\guillaume\\Data\\Mercator

is the temporary folder where the netcdf files from Mercator arrive every week and where they are converted into HDF5. Two HDF5 files are created per week. One for the Analyzed results of the past 7 days (labeled **ANA**). And one for the Predicted results of the coming 14 days (labeled **PRE**).

#### \\\\guillaume\\Aplica\\PreOp-Model\\GeneralData

is the folder where the interpolations from the Mercator grid and MM5 grids to the Mohid grids are made.

#### \\\\guillaume\\Aplica\\PreOp-Model\\WestIberia

is the folder where the post-processing is made to the results from the model every week.

#### \\\\192.168.23.151 or \\\\Dual\_64... \\Aplica\\PreOpModel

is the folder where the model is run every week.

#### ~\\\\kepler\\DataCenter\\Modelos\\Mercator-IST~ \\\\alexandria\\public\\modelos\\Mercator-IST

is the folder where the results in netcdf and in hdf5 from Mercator-Océan are stored for futur use.

#### \\\\data.mohid.com\\data

is the folder where the results in netcdf from MOHID are stored for publication via OpenDap.

## Scheduled jobs

Here's so far the approximate scheduling of the pending jobs on their respective locations:

### \\\\guillaume

1. 19:00 Thursday - Mercator Results mail reception
2. 20:00 Thursday - **Mercator-weekly**
3. 02:00 Friday - **Interpol-weekly**
4. 12:00 Friday - MM5 results
5. 17:05 Friday - **InterpolMM5-weekly**
6. 00:30 Monday - **PostRun-weekly**
7. 04:00 Monday - **PostRunHIDROMOD-weekly**

### \\\\192.168.23.151\\guillaume

```
0       2       *       *       6       /home/guillaume/Aplica/PreOp-Model/WestIberia_ET2_2K4/exe/CreateRun.sh
30      2       *       *       6       /home/guillaume/Aplica/PreOp-Model/WestIberia_ET2_2K4/exe/ MailRun.sh
0       2       *       *       1       /home/guillaume/Aplica/PreOp-Model/WestIberia_ET2_2K4/exe/MailRun.sh
0       4       *       *       *       /home/guillaume/Desktop/wgetwiki.sh
0       3       *       *       *       /home/guillaume/Aplica/PreOp-Model/cronjob.sh
```

## Actions StoryLine

### Mercator-Weekly

1. Calculate date and run ID
2. Access mercator ftp site and get all 63 netcdf files
3. Gunzip netcdf files
4. Gunzip last week last day netcdf file
5. convert to HDF5 the netcdf files for the analysis period
6. convert to HDF5 the netcdf files for the prediction period
7. Gzip back last week netcdf file
8. Gzip back the netcdf files
9. Transfer the netcdf files to **\\\\kepler**
10. Transfer the HDF5 files to **\\\\kepler**
11. Pass the run ID to the interpolation folder in **\\\\guillaume**
12. Send an acknowledgement email to Mercator

### Interpol-Weekly

1. Read the run ID
2. Glue the Analysis and Prediction HDF5 files to a Mercator HDF5 file
3. Interpolate the Mercator HDF5 to the Portugal grid producing a Portugal HDF5 file
4. Analyze the Portugal HDF5 and create a PortugalAnalyzed HDF5 file (barotropic flux)
5. Copy the Portugal HDF5 file to **\\\\192.168.23.151**
6. Copy the PortugalAnalyzed HDF5 file to **\\\\192.168.23.151**
7. Interpolate the Mercator HDF5 to the Estremadura grid producing a Estremadura HDF5 file
8. Analyze the Estremadura HDF5 and create a EstremaduraAnalyzed HDF5 file (barotropic flux) **NOT NECESSARY**
9. Copy the Estremadura HDF5 file to **\\\\192.168.23.151**
10. Copy the EstremaduraAnalyzed HDF5 file to **\\\\192.168.23.151**
11. Pass the run ID to the MM5 interpolation folder in **\\\\guillaume**
12. Send an acknowledgement email to administrator

### InterpolMM5-Weekly

1. Read the run ID
2. Glue this week MM5 HDF5 with last week MM5 file for domains D2 and D3
3. Patch D2 and D3 HDF5 files for the Portugal grid producing a Portugal HDF5 file
4. Copy the Portugal HDF5 file to **\\\\192.168.23.151**
5. Patch D2 and D3 HDF5 files for the Estremadura grid producing a Estremadura HDF5 file
6. Copy the Estremadura HDF5 file to **\\\\192.168.23.151**
7. Pass the run ID to the MM5 interpolation folder in **\\\\192.168.23.151**
8. Send an acknowledgement email to administrator

### CreateRun

1. Read the run ID
2. Calculate the START and END times
3. Create the dat files with the run ID
4. Run MohidWater as a forked process and returning output to a log file

### MailRun

1. Read the run ID
2. Tail the run log file to a status file
3. Send an email with the status file to the administrator

### PostRun-Weekly

1. Get the run ID file from **\\\\192.168.23.151**
2. Read the run ID
3. Get the log file from **\\\\192.168.23.151**
4. Match "*succefully*" string in log file
5. Copy the Portugal Hydrodynamic HDF5 file from **\\\\192.168.23.151**
6. Copy the Portugal WaterProperties HDF5 file from **\\\\192.168.23.151**
7. Convert to netcdf the Portugal Hydrodynamic HDF5 file
8. Copy the netcdf Portugal Hydrodynamic file to **\\\\data.mohid.com**
9. Convert to netcdf the Portugal WaterProperties HDF5 file
10. Copy the netcdf Portugal WaterProperties file to **\\\\data.mohid.com**
11. Copy the Estremadura Hydrodynamic HDF5 file from **\\\\192.168.23.151**
12. Copy the Estremadura WaterProperties HDF5 file from **\\\\192.168.23.151**
13. Convert to netcdf the Estremadura Hydrodynamic HDF5 file
14. Copy the netcdf Estremadura Hydrodynamic file to **\\\\data.mohid.com**
15. Convert to netcdf the Estremadura WaterProperties HDF5 file
16. Copy the netcdf Estremadura WaterProperties file to **\\\\data.mohid.com**
17. Send an acknowledgement email to administrator

### PostRunHIDROMOD-Weekly

1. Get the run ID file from **\\\\192.168.23.151**
2. Read the run ID
3. Get the log file from **\\\\192.168.23.151**
4. Match "*succefully*" string in log file
5. Copy the Portugal Hydrodynamic HDF5 file from **\\\\192.168.23.151** to **\\\\192.168.24.11\\c$\\www\\OperationalResults**
6. Copy the Portugal WaterProperties HDF5 file from **\\\\192.168.23.151** to **\\\\192.168.24.11\\c$\\www\\OperationalResults**
7. Copy the Portugal Atmosphere HDF5 file from **\\\\192.168.23.151** to **\\\\192.168.24.11\\c$\\www\\OperationalResults**
8. Send an acknowledgement email to administrator and end-users

## External References

- [data.mohid.com](http://data.mohid.com/index.xml),
- [model results in netcdf datasets](http://data.mohid.com/data.xml),
- [Access to a LAS](http://data.mohid.com/las.xml).
