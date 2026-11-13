# Mohid Statistics Analyser

The **Mohid Statistics Analyser** application calculates statistics and obtains graphs of comparison between two time series data.

Running options for this application are specified in an [interface](mohid-statistics-analyser.md).

## Introduction

For model validation or model comparison studies it is usefull to calculate comparison statistics and graphs. The Mohid Statistics Analyser allows the creation of such information for one or more properties with two site specific [time series](../../io/formats/time-series.md) files in a user-friendly way.

The properties which can be object of such analysis are the properties present in the time series. It is not required by the Mohid Statistics Analyser that the comparison is made for the same property in both time series, i.e., interproperties comparisons are possible. Per application run only one pair of time series can be analysed.

The information is obtained based on data located on a time window. That time window is defined as the time period for which data exist in both time series and any further constrain indicated by the user in the interface. If time instants in the two time series are different then the second time series data is linearly interpolated for the first time series data time instants.

The instant data in the time window is considered to obtain the all data information: "Statistics (All)" and "Graph (All)". This data is also preprocessed to derive the daily and monthly information: daily information ("Statistic (Daily)" and "Graph (Daily)") is obtained from the daily averages of the data and monthly information ("Statistic (Monthly)" and "Graph (Monthly)") is obtained from the monthly averages of the data.

The statistics which are calculated for each application run consist in:

- average of time series 1 data ("Observed Average");

- average of time series 2 data ("Modeled Average");

- bias;

- root mean square error (RMSE);

- Pearson product-moment correlation coefficient (R^2^);

- Nash-Sutcliffe coefficient of model efficiency.

The graphs are representations of the second time series data as a function of the first time series data - XY graphs. Several graph types and appearances can be selected.

**Typical use:**

Calculate statistics and graphics comparing the MOHID results with observed data.

**Data input requirements:**

Two [time series](../../io/formats/time-series.md) files containing data from different sources, generally one having observed data and another having model results.

**Ouput:**

Statistics and XY graphs, presenting the second time series data as a function of the first time series data, considering all data values for the time instants of the first time series, daily averages and monthly averages.

## User interface

The user interface is generally installed in the personal computer together with MOHID GIS and MOHID GUI and other applications, in the MOHID package.

The access to the interface executable file can be done through the Programs section of the Start Menu of the Desktop.

After access to the interface a **window** appears with **7 selectable tabs:**

- **Geral**: tab where the user makes the **input** of information required to run the Mohid Statistics Analyser application;

- **Statistic (All)**: tab presenting **output** statistics obtained considering all data values;

- **Statistic (Daily)**: tab presenting **output** statistics obtained from daily averages values;

- **Statistic (Monthly)**: tab presenting **output** statistics obtained from monthly averages values;

- **Graph (All)**: tab presenting the **output** XY graph of all data values;

- **Graph (Daily)**: tab presenting the **output** XY graph of daily averages values;

- **Graph (Monthly)**: tab presenting the **output** XY graph of monthly averages values.

In the next sections are presented the details about [Data input](mohid-statistics-analyser.md), [Statistics output](mohid-statistics-analyser.md) and [Graph output](mohid-statistics-analyser.md).

### Data input

The data input information is inserted by the user in the **Geral** tab, through the following steps:

- select the first time series, which will be considered the observed data and will define the time instants of both time series data, by writing the correspondent path in **Geral --> Observed Data --> Filename** or by clicking the respective **Browse** button and finding in the file in the computer folders; if the specified path does not exist then a window with an error message is generated;

- in the **Geral --> Observed Data** area appears the list of variables with data available in the selected time series; the user must then select the intended variable for the comparison;

- select the second time series, which will be considered the model data, by writing the correspondent path in **Geral --> Model Data --> Filename** or by clicking the respective **Browse** button and finding in the file in the computer folders; if the specified path does not exist then a window with an error message is generated;

- in the **Geral --> Model Data** area appears the list of variables with data available in the selected time series; the user must then select the intended variable for the comparison; this choice must be done carefully since the application does not require that the same variable must be selected in both time series;

- in the **Geral --> Parameter** area appears the time window where data is available for comparison; the user can restrict this time window to the desirable comparison interval using the available buttons; the time window can be reset at any time to the original one by clicking the **Parameter --> Reset Dates** button;

- click the **Geral --> Process...** button for input data processing.

During input data processing some messages for the user are generated and those appear together with a progress bar in the area bellow the **Process...** button. Error or warning messages can appear in separate windows during the processing, namely if the time window is badly specified (error message) or if a specific type of average data (daily or monthly) is not available.

The sucessfull processing terminates with a "Done!" message. In this case in the **Geral --> Show Time Series** presents for selection three buttons each causing the appearance of an individual window (**Time Serie Graph**) where the actual time series considered for comparison is presented in graph form:

- **All Values**: time series of values considered for Statistic (All) and Graph (All) tabs results;

- **Daily Values**: time series of values considered for Statistic (Daily) and Graph (Daily) tabs results (daily averages);

- **Monthly Values**: time series of values considered for Statistic (Monthly) and Graph (Monthly) tabs results (monthly averages).

In each time series graph window it is possible to change the features of the graph using the available buttons in the menu located above the graph, namely:

- **Chart Type**: opens a **Command and Options** window allowing the selection of several actions: **General**, **Border/Fill**, **Data Sheet**, **Type**, **Series Groups**, **Show/Hide**;

- **Undo**: undoes any graph change operation;

- **Show/Hide Legend**: selecting the button causes the appearance of the legend, clicking on it again makes the legend disapear;

- **Commands and Options**: opens a **Command and Options** window allowing the selection of several actions: **General**, **Border/Fill**, **Data Sheet**, **Type**, **Series Groups**, **Show/Hide**;

- **Chart Wizard**: opens a **Command and Options** window allowing the selection of actions concerning **Data Source**, **Data Sheet** and **Type**;

- **Help**: opens a help window where can be found information on how to change the graph features using the buttons referred before.

Four additional buttons, located beneath the graph, command the following actions:

- **CheckBox1**: checking it allows the selection in **Window Width 

$$
days
$$

** the number of days to appear in the graph and in the ruler the location in the time interval of that number of days;

- **Save XML**: saves the graph in XML code in the path provided by the user;

- **Save Image**: saves the graph as a GIF or JPEG file in the path provided by the user;

- **Exit**: causes the closing of the **Time Serie Graph** window.

Back to the Mohid Statistics Analyser Geral tab, the lower-right button **Exit** exits the Mohid Statistics Analyser window.

### Statistics output

The Mohid Statistics Analyser statistics output is organized in the **Statistic (All)**, **Statistic (Daily)** and **Statistic (Monthly)** tabs, for statistics considering all data, daily averages and monthly averages, respectively.

This output is only produced after the data input by the user in the **Geral** tab and the successful data processing.

In each of the statistics tabs the information content is similar:

- **Observed Average** (average of the first time series);

- **Modeled Average** (average of the second time series);

- **Bias**;

- **RMSE** (root mean square error);

- **R^2^** (Pearson product-moment correlation coefficient);

- **Model Efficiency** (Nash-Sutcliffe coefficient).

For each of these statistic parameters is presented the calculated value, in the grey area, and a small description of the meaning and the interpretation of the parameter, in the blue area.

Beneath the parameters information are presented some notes on the approach used to calculate the statistics.

The lower-right button **Exit** in each Statistics tab exits the Mohid Statistics Analyser window.

### Graphic output

The Mohid Statistics Analyser graphic output is organized in the **Graph (All)**, **Graph (Daily)** and **Graph (Monthly)** tabs, for graphs considering all data, daily averages and monthly averages, respectively.

This output is only produced after the data input by the user in the **Geral** tab and the successful data processing.

In each of the graphs tabs the information content is similar: a XY point graph of second time series data (**Model Data**) as a function of the first time series data (**Observed Data**).

It is possible to change the features of the graph using the available buttons in the menu located above the graph, namely:

- **Chart Type**: opens a **Command and Options** window allowing the selection of several actions: **General**, **Border/Fill**, **Data Sheet**, **Type**, **Series Groups**, **Show/Hide**;

- **Undo**: undoes any graph change operation;

- **Show/Hide Legend**: selecting the button causes the appearance of the legend, clicking on it again makes the legend disapear;

- **Commands and Options**: opens a **Command and Options** window allowing the selection of several actions: **General**, **Border/Fill**, **Data Sheet**, **Type**, **Series Groups**, **Show/Hide**;

- **Chart Wizard**: opens a **Command and Options** window allowing the selection of actions concerning **Data Source**, **Data Sheet** and **Type**;

- **Help**: opens a help window where can be found information on how to change the graph features using the buttons referred before.

The lower-right button **Exit** in each Graph tab exits the Mohid Statistics Analyser window.

The Mohid Statistics Analyser application does not allow the saving of the Graph tabs images.
