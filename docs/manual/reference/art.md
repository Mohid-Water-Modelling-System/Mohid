# ART

## About

ART (Automatic Running Tool) is a tool for managing operational models, this tool is programmed in Python and adapted from a previous C++ version. The new python version allows to operate models in windows and linux environemnts. The code can be found at XXXXXXXXXX. A short manual and the description of the keywords can be found below.

## Software needed

##### YAML Keywords

#### ARTCONFIG

- MAIN\_PATH - MANDATORY - path to the main project, needs to exist.
- OPERATIONAL\_MODE - if set to \`1\` it sets \`ART\` mode to run as an operational model, it uses the current day on the machine and \`REF\_DAYS\_TO\_START\` to calculate the \`START\_DATE\` and uses \`DAYS\_PER\_RUN\` and \`NUMBER\_OF\_RUNS\` to calculate \`END\_DATE\`. If set to \`0\` or non existent it defaults to a single run (using \`START\_DATE\` and \`END\_DATE\` parameters instead of calculating them).
- RUN\_PREPROCESSIN - if set to \`1\` it will run all programs defined in the \`PRE\_PROCESSING\` block. If set to \`0\` or non existent it will default to not running them.
- DAYS\_PER\_RUN
- REF\_DAYS\_TO\_START - reference day relative to current day on the machine. It sets the \`START\_DATE\` for \`OPERATIONAL\_MODE\`.
- NUMBER\_OF\_RUNS
- MODULE - defines which program \`ART\` is going to use to calculate the model. Allowed values: \`mohid\`, \`ww3\`, \`wrf\`.
- RUN\_SIMULATION\_MODULE - if set to \`1\` it will run the program specified in the \`MODULE\` parameter. If set to \`0\` or non existent it will default to not running it.
- RUN\_SIMULATION\_POST - if set to \`1\` it will run all programs defined in the \`POST\_PROCESSING\` block. If set to \`0\` or non existent it will default to not running them.
- START\_DATE - start date to calculate model in a non-operational run, ignored if \`OPERATIONAL\_MODE\` is set.
- END\_DATE - end date to calculate model in a non-operational run, ignored if \`OPERATIONAL\_MODE\` is set.
- OUTPUT\_TO\_FILE - if set to \`1\` it will redirect the output of the console (\`ART\`) to a file specified by \`OUTPUT\_FILEPATH\` parameter. If \`0\` or non existent it will only output to the console where it is running.
- OUTPUT\_FILEPATH - specifies the file path for which \`ART\` will output. If non existent and \`OUTPUT\_FILE\` is set to \`1\` it will default to the path where \`ART\` is running.
- SEND\_EMAIL - if set to 1 it will enable sending emails upon the completion of \`MOHID\` runs. If 0 it disables that option.
- EMAIL - list of strings that specify the address(es) of the person the email is sent to. If \`SEND\_EMAIL\` is set to 1, it has to exist, and it has to be non-empty.

#### MOHID

- MAX\_TIME - maximum time allowed for \`MOHID\` to run, if non existent it will default to \`40000\`.
- OUTPUT\_TO\_FILE - if set to \`1\` it will redirect the output of \`MOHID\` to a file specified by \`OUTPUT\_FILEPATH\` parameter. If \`0\` or non existent it will only output to the console where it is running.
- OUTPUT\_FILEPATH - specifies the file path for which \`MOHID\` will output. If non existent and \`OUTPUT\_FILE\` is set to \`1\` it will default to the path where \`MOHID\` is running.
- EXE\_PATH - MANDATORY - specifies the path to where the \`MOHID\` execution file is
- MPI
- ENABLE - if set to \`1\` it will enable \`MPI\` on current run. If \`0\` or non existent it will default to ignore all the \`MPI\` block
- EXE\_PATH - MANDATORY IF \`ENABLE\`- specifies the path to where the MOHID execution file is.
- TOTAL\_PROCESSORS - number of total processors cores that \`MPI\` will use

#### model\_name

```
  - NAME
  - PATH - path relative to `MAIN_PATH`
  - GRID_FILE
  - DT
  - STORAGE_PATH
  - RESULTS_LIST - defined as a list `['Hydrodynamic_surface', 'WaterProperties']`. When used makes so that the Backup of the results is only made if the file name is in the list. If this parameter is not present the behavior will default to backup all files.
  - HAS_SOLUTION_FROM_FILE
  - DISCHARGES
     - discharge block name
        - ENABLE
        - PATH
        - DATE_FORMAT
  - OBC
     - obc block name
        - ENABLE
        - FILE_TYPE - file type of the OBC file (usually hdf5 or netcdf). When not defined it defaults to `hdf5`.
        - SIMULATED_DAYS
        - SUB_FOLDERS - if set to `1` it searches it assumes that `OBC` `WORKPATH` is structured with `SUB_FOLDERS` for `year` if set to `2` it has `SUB_FOLDERS` `month` and `year` and if set to `3` it has `year` , `month`, and `day`. If `0` or non defined it defaults to assuming no structure inside `OBC` `WORKPATH`.
        - DATE_IN_FILENAME -  if set to `1` assumes date is in the files' names (uses `dateFormat`). If `0` or non defined it defaults to not having date on the obc files.
        - DATE_FORMAT - determines the date format for the obc files' names. If not defined it defaults to `YYYYMMDD`
        - FILES - list of files you want from the OBC workpath `['Hydrodynamic', 'WaterProperties']` it can be costumized with date variables `%Yi, %Mi, %di` for initial year, initial month and initial day and also  `%Yf, %Mf, %df` for final year, final month and final day.
        Example:
        `['Hydrodynamic_%Yi-%Mi-%di']` will become `Hydrodynamic_2019-09-12`.
        - WORKPATH
        - METEO
           - ENABLE
           - MODELS
              - uniqueId
                 - NAME
                 - SIMULATED_DAYS
                 - FILENAME_FROM_MODEL
                 - WORKPATH
                 - DATE_FORMAT
                 - FILE_TYPE
```

#### MODEL.DAT

```
   Arguments used here will be inserted or changed in the model.dat of its domain.
   - MAXDT
   - GMTREFERENCE
   - DT_PREDICTION_INTERVAL
```

### PREPROCESSING

```
  - name of block - must be unique within the `PREPROCESSING` block
     - RUN
     - WORKING_DIRECTORY - sets the `WORKING_DIRECTORY` of the pre_processing tool to the one specified.
     - DAT_DATE_CHANGE - if the .dat of the script needs to change its START and END set this parameter to `1` and define `CONFIG_FILEPATH`.
     - CONFIG_FILEPATH - path to where the .dat file is. Only necessary when `DAT_DATE_CHANGE` is enabled.
     - EXE_PATH - mandatory - path to the executable that you want to run before the simulation.
     - FLAGS - any flags or arguments you must give in the command line to the executable
       ex. `ls -l -a /samba/lusitania` → flags: "-l -a /samba/lusitania"
     - OUTPUT_TO_FILE
     - OUTPUT_FILE_PATH
```

### POSTPROCESSING

```
  - name of block - must be unique within the preprocessing block
     - RUN
     - WORKING_DIRECTORY - sets the `WORKING_DIRECTORY` of the post_processing  tool to the one specified.
     - DAT_DATE_CHANGE - if the .dat of the script needs to change it START and END set this parameter to `1` and define `CONFIG_FILEPATH`.
     - CONFIG_FILEPATH - path to where the .dat file is. Only necessary when `DAT_DATE_CHANGE` is enabled.
     - EXE_PATH - mandatory - path to the executable that you want to run after the simutaltion. Can also be used to run a shell command such as `ls`. If you want to run a script such as `flags` parameter.
      Example: `exePath` : `python3` ; `flags` :`/pathToScript/script.py`
     - FLAGS - any flags or arguments you must give in the command line to the executable
     ex. `ls -l -a /samba/lusitania` → flags: "-l -a /samba/lusitania"
     - OUTPUT_TO_FILE
     - OUTPUT_FILEPATH
```
