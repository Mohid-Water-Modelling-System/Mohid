# Calibration/Validation in Mohid Land

## Resume

This page was written intended to give the main insights to calibration/validation in Mohid Land for hydrology. Since the model is physically based the main guidances are linked to giving the model the right precipitation, evapotranspiration, geometry in sections or impermeability matrix representing the area available for infiltration. This section describes the steps for conducting a procedure for calibration and should be done sequencially since the first steps exist to assure that the last are fine tunings and not deviations needed to compensate the lack of work in the first steps.

There are 4 main steps:

1 - Validate rain data since is the main driver for all the water cycle (long term - annual and monthly - and short term - during floods). Also validate other meteorology data (temperature, relative humidity, radiation, wind velocity) important for reference evapotranspiration.

2 - Validate actual evapotranspiration or river flow in long term since evapotranspiration is the link between rain and river flow.

3 - Validate floods rise and fall (short term dynamics) that is composed mainly by surface water.

4 - Validate groundwater flow that has it most component after the flood (after aquifer rise and transport trough the porous media)

## Before starting

Model lack of accurancy usually is simpton of three things:

- or model does not describe fully the system physics/relations

- or model data provided (initial or boundary conditions) is not good enough to describe its state or variability (spatially or temporally)

- or both,

If we agree that Mohid Land model has the main physical processes that describe watershed hydrodynamics than one should have in mind that model results accurancy will suffer from the given data description of reality. And those data consist of:

- Digital Terrain model that will affect land and river slopes, drainage direction and accumulation, etc

- Land Use map that affects the type of vegetation growing, runoff resistance to flow, impermeabilization, etc

- Soil Map and van Genuchten parameters that affect porous media flow, infiltration, etc

- Rain data that is the main driver of the water cycle and can have high spatial and temporal variation and usually prone to errors (collection errors, publication errors, data manipulation errors) that may impact severely the results.

In the next chapters it will be addressed the calibration/validation method that should be used, assuming that in terms of DTM, soil and land use map the user is using the best available data.

## Step 1 - Validate rain and meteorology

Since rain is the main driver of the water cycle, data gaps or data errors will compromise results. For example if looking to simulate floods there can not be a flood in the model if it does not rain in that period, or model will generate really high floods if abnormally high rain is inputed.

It is very important that the rain used for feeding the model is the best one. And the best is the rain that represents the rain regime and events occured in the watershed. With that said it is important that if the watershed has areas with higher rain intensity or lower that they are represented by rain stations so that the difference will be taken in account by the model.

So **What stations to use and how to trust in data?**

The first things to do is to look at rain stations. Collect all the data in the watershed but also in neighbour watersheds and do annual and monthly accumulations of data with the number of sampling days per month and per year. When graphing the accumulated annual and monthly precipitations, the stations (or periods) with bad data will pop up trough low sampling days or curves in graphs that are not correlated to other stations.

In annual accumulation data, stations should be correlated and the highest or lowest rain stations should maintain the ranking. In annual accumulated rain graphs also is possible to see the areas with higher rain intensities and lower and the number of stations that need to use in order to capture the spectrum of precipitation regimes. It is also a good help to use maps of annual precipitations to check the consistency of the data or the lack of a precipitation station for a high or low precipitation regime; in the case of lacking stations inside the watershed to represent the precipitation regimes it should be taken from a station outside the watershed.

So removing stations with uncorrelated or inconsistent data and using the stations that are consistent between each other and represent the precipitation regimes is the method right to select precipitation data.

In the case that looking for fast floods than additionaly the precipitation stations should have sub-daily data and rain events should occur at the same period of the recorded flow or level, if not the model obviously will not be able to represent the floods.

For other variables as radiation, relative humidity, temperature, wind velocity that interfere for instance on reference evapotranspiration computation, one should also compare its data to neighbour stations in annual and monthly analysis to be sure that they are consistent. Since this meteorology is less spatial variable than rain and has a more regional variation, the stations values should be very similar between stations. Stations or periods not consistent should be eliminated and because of the regional variation data from stations more distant from watershed may be considered.

## Step 2 - Validate evapotranspiration or river flow in long term simulation

After the rain is consistent and representative of the precipitation regime(s) in the watershed it is needed to validate the long term water fluxes and the weight of evapotanspiration.

If actual evapotranspiration data is available (e.g. EO data) than the comparison should be done directly to that data.

Wether or not actual evapotranspiration data is available it should always be compared model flow to collected data (start with yearly and monthly analysis). If rain is consistent and flow in river has a good agreement in monthly and annual annalysis than the actual evapotranspiration is being correctly estimated in the long term.

**Things to check:**

Reference Evapotranspiration is computed from atmosphere properties and should be checked with other sources of data.

**Parameters:**

Kc - relates reference evapotranspiration (from alfafa) to crop transpiration (specific of crop and plant density)

Feddes head limits - represent the crop stress to water contents

## Step 3 - Validate floods (surface water)

Important: Remind that if changing this parameters it may change evapotranspiration and river flow and after this step the long term validation discussed in the previous section should be redone.

After rain check (in long term and flood periods) in first step and evapotranspiration/river flow check for long term (annual, monthly) than model is predicting the right amount of evapotranspiration and river flow for monthly analysys. However for floods that may even rise and fall in one day additional calibration may be needed.

As a rule of thumb impermeabilization is more important for controlling first floods of the year or a first rise of several floods in a row and evapotranspiration (Kc or Feddes parameters) is more important after first floods or first rise. This is explained by the fact that first floods occur with soil not dry so it has ability to infiltrate water and evapotranspiration is not so important since excess saturation does not occur and the way to control the water that arrives to the river is trough impermeabilization. On the other side, for second floods soil is already saturated and more or less impermeabilization is the same since the soil is already behaving as impermeable and the most important factor to control water that goes to river is evapotranspiration that if higher removes more water from soil and leaves space for imfiltration and evapotranspiration, water that will not arrive to the river.

### Flood Rise

#### Model is delayed

**Description:**

If the model is not getting a flood rise but data has, check that rain exists at that period. If rain exists check that the water content in soil increased. What happened is that in the model water infiltrated and did not created the flood rise at that point and only after filling soil the floods generated (delayed to data). This may happen usually at first peaks were soil is not dry.

**Solution:**

1 - Check if the rain itensity is right (with other stations) because runoff may be generated not only by saturation excess (filling soil) but also by high rain itensity (infiltration excess)

2 - Usually the first floods because soil is dry (or drier) do not occur from saturation excess or infiltration excess (if not really high intensities) but impermeabilization is the main control factor. So if the watershed is urban it should be described the impermeable fraction and the storm water structures; if the watershed is rural than impermeability may occur beacause of soil sealing.

Also the drainage network may have low resolution (only main rivers are described) and in reality concentrated flow (flow that occurs concentrated in a width lower than cell size) may occurr. So, before anything increase your drainage network resolution and see if the first peaks start to appear. If yes it means that in reality exists concetrated flow inside the cells that can only be simulates with drianage network (by the way that was the reason why drainage network module was created, because if the DTM and cell resolution were of the order of 1m or less there was no need for it).

**Parameters:**

- soil impermeability - a grid or HDF with impermeability ratio

- storm water infiltration - infiltration velocity to storm water drainage system

- storm water velocity - water velocity in storm water drainage system

#### Model is earlier

**Description:**

If the model is getting a flood rise faster than data.

**Solution:**

1 - Check if the rain is right (with other stations) and is not previous to the flow data.

2 - Check if impermeabilization or storm water are not creating a flood to soon (to much impermeabilization or to fast storm water velocity).

**Parameters:**

- soil impermeability - a grid or HDF with impermeability ratio

- storm water infiltration - infiltration velocity to storm water drainage system

- storm water velocity - water velocity in storm water drainage system

### Floods After Soil Saturation

As stated above second floods or floods after soil is saturated are less dependent on impermeabilization factor (soil is behavinf already as impermeable) and evapotranspiration is the controling factor since it may increase or decrease soil water content and its ability to infiltrate more or less water that will not get into the river.

**Parameters:**

- Kc

- Feddes Heads

### Flood Fall

#### Model has a faster decay

**Description:**

If model fall is faster than data

**Solution:**

1 - Check if changing the channel manning (higher mannings upstream) because of bottom type may be enough to delay water that is upstream.

2 - Check if the section is right with the measuring local. Check the section geometry including slope that is an important factor to remove water fast if it is bigger that it should. Also check for section interruptions (or sedimentation) as sections may become smaller in particular areas and the water may only flows trough a part of the section. Sections can be changed in Mohid Land by hand.

**Parameters:**

- Channel manning

- Section geometry (including bottom slope that is not a parameter but computed from adjacent nodes bottom levels).

#### Model has a slower decay

**Description:**

If model fall is slower than data

**Solution:**

1 - Check if lowering the channel manning reduces the effect because it may be creating a big friction and water level rises, taking more time for the flood.

2 - Check if the section is right with the measuring local. Check the section geometry including slope that is an important factor to remove water fast if it is slower that it should. Also check for section changes as widening or deepening. Sections can be changed in Mohid Land by hand.

**Parameters:**

- Channel manning

- Section geometry (including bottom slope that is not a parameter but computed from adjacent nodes bottom levels).

## Step 4 - Validate baseflow

Important: Remind that if changing this parameters it may change evapotranspiration and river flow and after this step the long term validation discussed in the section 3 should be redone.

After the flood rise and fall (surface water) the most of the water that arrives in river is from groundwater, because aquifer rise as response to rain but the velocity in porous media is much slower so the water arrives to the river after the surface water.

### Model Baseflow is low and falls fast

**Description:**

The model results hydrogram falls fast after the flood and almost no base flow is maintained.

**Solution:**

The groundwater flow needs to be enhanced. Mohid Land computes groundwater flow based on the Hydraulic gradient between river level and aquifer level and the aquifer level is computed up to where saturation occurs. Model has cells and fluxes between cells (specially in periods of high flows) may lead to a situation of soil very near saturation but not completly saturated and so the aquifer level is artificailly lowered so as the flux between river and soil.

Give a water content range - from saturation - to check for groundwater level. The default value is 1E-15 difference to saturation water content that creates a high aquifer level fluctuation (an mainly reduction).

**Parameters:**

- groundwater level range for saturation - range from saturation point when searching for aquifer level.

- horizontal conductivity factor - horizontal conductivity multiplying factor (multiplying for cell computed conductivity).

## Links

[Mohid Land](../land/mohid-land.md)
