# Mohid Time Series Editor

## Overview

Mohid Time Series Editor is a graphical user interface which allows the user to visualize in a quick way [time series](../../io/formats/time-series.md) data required or produced by the Mohid numerical programs.

Usually Mohid Time Series Editor will be launched directly from Mohid GUI, opening directly the file chosen by the user (see [Switching to Post Processing Mode](mohid-gui.md)).

Mohid Time Series Editor is written in VB.NET and uses Office Web Components to display the results.

## Open a file

By selecting **"File -> Open"** from MOHID Time Series Editor’s main menu a new time series is added to the list view on the left side of MOHID Time Series Editor’s main window.

## Selecting series to display

After opening a new file a dialog window appears. Here the user has the option to select which time series are to be displayed. The process of opening files and selecting series must be repeated for all time series to display.

## Display data

By pressing the **"Draw Chart"** button at the lower right corner in MOHID Time Series Editor’s main window, the Chart Window which displays the selected time series will appear. By selecting commands from the toolbar at the top of the Chart Window, the chart can be manipulated in a similar way as in MS Excel. By selecting the **"Save Image"** button the current image is saved.

## Save entire graphic

The user is allowed to save an entire graphic (values and settings) by writing it to XML format. To do this press the SaveXML and select the file name. The entire graphic can be loaded later by selecting **"File -> Open XML Time Serie"**.
