# David ToDoList

## Article Tasks



$$
Task
$$

 

$$
Priority
$$

 

$$
MonthToFinish
$$



1. Validate Model with analitical solutions - porous media and runoff
2. Validate Model with erosion/deposition data (Jose)
3. Compare Vegetation growth and erosion with SWAT in subbasin of Save and Enxoe
4. Colocar modelo pesticidas MOHID
5. Article 1 - part 1 (David Brito) Hidrology of MOHID Land and validation - analitic; part 2 (David Bailly) - Erosion and validation (José data)
6. Article 2 - SWAT and Mohid Land
7. Article 3 - application in Save (D Bailly)
8. Article 4 - application in Enxoé (D Brito) - functioning of basin in one year and in events
9. Article 5 - pesticides

## Coding Tasks

### General

1. Concept a structure for using basin as the model that computes all interface fluxes and gives them to the modules that need them - fluxes are computed at a given state after or before all internal module fluxes. Think as a property flux (what concentration use).
2. Concept a module to handle properties? for instance one property defined for all the modules.

### Geometry

1. Rewrite of bottom file should be a warning and should be used the good file (and for instance saved the user original \_.old in order not to create inconsistency between Porous Media data file (BOTTOM\_FILE) and file effectively used).
2. When geometry depth is lower than maximum domain depth the warning should be more explicit to say to the user to increase geometry depth.

### Module Basin

1. Compute Mass Balance for each property CHECK
2. Take WATER\_COLUMN\_COEF from module and data file CHECK
3. Water and Mass from rain should be updated after atmosphere processes and so it does not need to compute infilcolconc for porous media properties (it is the water column) CHECK
4. If evaporation removes all leaf water then mass should be kept - use mass instead of concentration - CHECK
5. If evaporation removes all water column then mass should be kept - use mass instead of concentration.

### Module Porous Media

1. Soil ID's in <beginsoiltype> <endsoiltype> should not need to beconsecutive.
2. Remove warnings from screen (river section lower than soil) and write them for example in error and messages. Left just one comment that something happended
3. Chage connection with river for layer based CHECK
4. Infiltration - calibrate runoff
5. Convergence criteria by % of domain cells instead of one cell
6. If convergence criteria is not met - dt is cutted and solution should go back to last known good configuration instead of going back to beggining of dt
7. Distribute oversaturation differently - use oversaturation to Head field - physic approach
8. Water that exists the soil in sub VerticalContinuity in a river cell should go to river directly and not to runoff first.

### Module PorousMediaProperties

1. Test interaction with river CHECK
2. Construct log to screen with properties CHECK
3. Get transpiration flux in to advection diffusion CHECK
4. Nutrient uptake may not be linked to flow\*conc and in porous media properties it is uptaked as sink. But the time that vegetation checkedthe amount and the time is taken may pass one day.
5. Restructure routine VegetationInterfaceFluxes to make it faster and straightforward
6. Implement SedimentQuality and test it - CHECK
7. Compare solute transport to numeric solution

### Module RunoffProperties

1. Construct bottom fluxes for properties - follow "BOTTOM\_FLUXES" in DrainageNetwork. Property will have a water matrix and bottom matrix. First approach can be only transport, concentrate when part of water infiltrates and deposition when all water infiltrates. Then after, erosion and deposition can be programmed (see Drainage Network again). CHECK
2. Construct partition between dissolved and particulate - follow "PARTICULATE" in water properties. CHECK. Need tests
3. Construct log to screen with properties CHECK
4. Add the Possibility to use implicit advection diffusion (thomas3D). CHECK
5. Get discharges flux and boundary flux in to advection diffusion
6. Create rills for erosion/deposition (rill cells and area for flux defined with drainage accumulation computation)

### Module Runoff

1. Construct, modify and kill timeseries CHECK

### Module Vegetation

1. Remove warnings from screen (root depth lower than soil) and write them for example in error and messages. Left just one comment that something happended
2. In the simplest case where there is no gowth model there is no need for <begin\_agricultural\_practices>....
3. Solve Continuous computation problem in properties - CHECK
4. Nutrient uptake may not be linked to flow\*conc (SWAT based) and in porous media properties it is uptaked as sink. CHECK - removed it from advection diffusion routine and put it as sink before transport computations.
5. the time that vegetation checked nutrient amount and the time is taken may pass one day and transport may remove the material. Think about that. CHECK for now if concentration negative mass creation matrix produced and warning print.
6. Code Fertilization events as pesticide events (for now only autofertilization) - CHECK
7. Code Rotation with hdf5 - CHECK
8. Code a secondary plant type that grows after main plant death - may be unnecessary with rotation. think about that.
9. Code the excess water stress in growth model
10. Create DB for vegetation and automatically generate data whitout writting in vegetation data file. Create a file to check what was read. Use David Bailly Work to write in on file - CHECK
11. Give to porous media properties fluxes in g/s
12. How to account leaf flux to soil - for dynamic manning roughness?
13. Check equation of SWAT to separate potential evaporation using LAI to see if is the same as used in Basin
14. Uptake also ammonia?

### Module Drainage network

1. Burn in drainage network - from shape or line
2. Code evapotranspiration in river to remove excess water present. when water depth reduces in summer and plants overtake river. Almost CHECK.
3. Unit test with 2 delineations - 1 outlet at the midle and next delineate the basin there and compare with the new outlet

## Project Tasks

Excel file
