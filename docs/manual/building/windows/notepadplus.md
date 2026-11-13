# Notepadplus

To use syntax highlighting and -folding in the editor notepad++ follow these steps:

1.) Download notepad++ ([http://notepad-plus.sourceforge.net/](http://notepad-plus.sourceforge.net/) you can choose the language there)

2.) install it

3.) copy the file userDefineLang.xml given below to the application data folder for notepad++ e.g. C:\\Dokumente und Einstellungen\\Username\\Anwendungsdaten\\Notepad++

4.) Start MohidGui (open a project)

5.) go to /tools/options

6.) choose "Data Files/ Tools"

7.) specify the path to notepad++ (e.g. C:\\Programme\\Notepad++\\notepad++.exe)

## userdefineLang.xml

```
<NotepadPlus>
   <UserLang name="Mohid.dat" ext="dat xyz">
       <Settings>
           <Global caseIgnored="no" />
           <TreatAsSymbol comment="no" commentLine="yes" />
           <Prefix words1="no" words2="no" words3="no" words4="no" />
       </Settings>
       <KeywordLists>
           <Keywords name="Delimiters">000000</Keywords>
           <Keywords name="Folder+">begindischarge BeginOrigin BeginProperty beginproperty BeginYY BeginXX BeginGridData2D BeginTimeSerie BeginNodeTimeSerie begin_xyz BeginBathymetry begin_rugosity begin_erosion_rate begindomain beginvegetationtype begin_property beginwaterlevel beginhorizon beginsoiltype begin_station</Keywords>
           <Keywords name="Folder-">enddischarge EndOrigin EndProperty endproperty EndYY endwaterlevel endhorizon endsoiltype end_station EndXX EndGridData2D EndTimeSerie EndNodeTimeSerie end_xyz EndBathymetry end_rugosity end_erosion_rate enddomain endvegetationtype end_property</Keywords>
           <Keywords name="Operators">- . : < ></Keywords>
           <Keywords name="Comment">1 2 0!</Keywords>
           <Keywords name="Words1"></Keywords>
           <Keywords name="Words2"></Keywords>
           <Keywords name="Words3">NAME UNITS DEFAULTVALUE DESCRIPTION NODE_ID DATA_BASE_FILE FLOW_COLUMN DT_OUTPUT_TIME MAX_BUFFER_SIZE Based on Drainage Network file Based on XYZ file NETWORK_FILE CHECK_NODES CHECK_REACHES HYDRODYNAMIC_APROX NUMERICAL_SCHEME MASS_ERR GLOBAL_MANNING MIN_WATER_DEPTH INITIAL_WATER_DEPTH XS_CALC DOWNSTREAM_BOUNDARY DISCHARGES STABILIZE MAX_ITERATIONS DT_FACTOR MAX_DT_FLOOD CONTINUOUS TIME_SERIE_BY_NODES TIME_SERIE_LOCATION TRANSMISSION_LOSSES HYDRAULIC_CONDUCTIVITY OUTPUT_NAME INITIALIZATION_METHOD DEFAULT_VALUE MIN_VALUE WATER_QUALITY TIME_SERIE BOTTOM_FLUXES EROSION CRIT_SS_EROSION EROSION_COEF DEPOSITION CRIT_SS_DEPOSITION CHS WS_TYPE WS_VALUE KL KL1 ML M START END DT VARIABLEDT GMTREFERENCE SPLITTING SOLARRADIATION_FILE SOLARRADIATION_COLUMN AIR_TEMPERATURE_FILE AIR_TEMPERATURE_COLUMN CLOUD_COVER_FILE CLOUD_COVER_COLUMN RELATIVE_HUMIDITY_FILE RELATIVE_HUMIDITY_COLUMN WIND_SPEED_FILE WIND_SPEED_COLUMN EVOLUTION CORIOLIS VOLUMEVARIATION HORIZONTALADVECTION HORIZONTALDIFFUSION INITIAL_ELEVATION INITIAL_VEL_U INITIAL_VEL_V VERTICALADVECTION VERTICALDIFFUSION VELTANGENTIALBOUNDARY VELNORMALBOUNDARY NULL_BOUND_HORADV TIDEPOTENTIAL WIND ATM_PRESSURE SURFACEWATERFLUX HMIN_CHEZY VMIN_CHEZY DATA_ASSIMILATION IMPLICIT_HORADVECTION RECORDING MOMENTUM_DISCHARGE SLOWSTART HMIN_ADVECTION LOCAL_DENSITY CONSERVATIVE_HOR_DIF BIHARMONIC BIHARMONIC_COEF SUBMODEL MISSING_NULL DEADZONE BAROCLINIC_RADIATION LOCAL_SOLUTION BOTTOMWATERFLUX SLIPPING_CONDITION NONHYDROSTATIC STATISTICS TIDE ENTERING_WAVE RADIATION WATER_DISCHARGES RESIDUAL ENERGY UP_CENTER IMPLICIT_VERTADVECTION IMPLICIT_VERTDIFFUSION DISCRETIZATION UPSTREAM BAROCLINIC BRFORCE OUTPUT_TIME I_CELL J_CELL DEPTHBYGRID K_CELL CONSTANT_CONC TIME_SERIE_COLUMN REFERENCE_DENSITY OUTPUT_HDF BOUNDARY_CONDITION IS_COEF PARTICULATE BOUNDARY_INITIALIZATION DEFAULTBOUNDARY ADVECTION_DIFFUSION ADVECTION_H_IMP_EXP VERTICAL_MOVEMENT SURFACE_FLUXES PARTITION WATER_FLUXES MASS_LIMITATION MASS_MIN BOX_TIME_SERIE OLD VISCOSITY_V VISCOSITY_H SHADING_FACTOR REMOVE_OVERTOP AERATION_METHOD FRACTION_SEDIMENT	 MINIMUM_SLOPE COMPUTE_LOAD	 SUMTOTALCONC LOCALIZATION_I LOCALIZATION_J LOCALIZATION_K</Keywords>
           <Keywords name="Words4"></Keywords>
       </KeywordLists>
       <Styles>
           <WordsStyle name="DEFAULT" styleID="11" fgColor="000080" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="FOLDEROPEN" styleID="12" fgColor="008000" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="FOLDERCLOSE" styleID="13" fgColor="008040" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="KEYWORD1" styleID="5" fgColor="0000FF" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="KEYWORD2" styleID="6" fgColor="008080" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="KEYWORD3" styleID="7" fgColor="004080" bgColor="FFFFFF" fontName="" fontStyle="1" fontSize="10" />
           <WordsStyle name="KEYWORD4" styleID="8" fgColor="FF8040" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT" styleID="1" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT LINE" styleID="2" fgColor="C0C0C0" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="NUMBER" styleID="4" fgColor="FF8080" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="OPERATOR" styleID="10" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER1" styleID="14" fgColor="FF0000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER2" styleID="15" fgColor="FFFF80" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER3" styleID="16" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
       </Styles>
   </UserLang>
   <UserLang name="Mohid.dnt" ext="dnt">
       <Settings>
           <Global caseIgnored="no" />
           <TreatAsSymbol comment="no" commentLine="no" />
           <Prefix words1="no" words2="no" words3="no" words4="no" />
       </Settings>
       <KeywordLists>
           <Keywords name="Delimiters">000000</Keywords>
           <Keywords name="Folder+"><BeginNode> <BeginReach></Keywords>
           <Keywords name="Folder-"><EndNode> <EndReach></Keywords>
           <Keywords name="Operators">.</Keywords>
           <Keywords name="Comment"></Keywords>
           <Keywords name="Words1">BeginNode EndNode BeginReach EndReach</Keywords>
           <Keywords name="Words2">ID GRID_I GRID_J COORDINATES CROSS_SECTION_TYPE DRAINED_AREA BOTTOM_LEVEL BOTTOM_WIDTH TOP_WIDTH HEIGHT</Keywords>
           <Keywords name="Words3">DOWNSTREAM_NODE UPSTREAM_NODE</Keywords>
           <Keywords name="Words4"></Keywords>
       </KeywordLists>
       <Styles>
           <WordsStyle name="DEFAULT" styleID="11" fgColor="004000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="FOLDEROPEN" styleID="12" fgColor="008080" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="FOLDERCLOSE" styleID="13" fgColor="008080" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="KEYWORD1" styleID="5" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="KEYWORD2" styleID="6" fgColor="808080" bgColor="FFFFFF" fontName="" fontStyle="2" />
           <WordsStyle name="KEYWORD3" styleID="7" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="4" />
           <WordsStyle name="KEYWORD4" styleID="8" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT" styleID="1" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT LINE" styleID="2" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="NUMBER" styleID="4" fgColor="FF80C0" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="OPERATOR" styleID="10" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER1" styleID="14" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER2" styleID="15" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER3" styleID="16" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
       </Styles>
   </UserLang>
   <UserLang name="Mohid.mgp" ext="mgp">
       <Settings>
           <Global caseIgnored="no" />
           <TreatAsSymbol comment="no" commentLine="no" />
           <Prefix words1="no" words2="no" words3="no" words4="no" />
       </Settings>
       <KeywordLists>
           <Keywords name="Delimiters">000000</Keywords>
           <Keywords name="Folder+">begin_view begin_bottom_legend begin_bottom_legend_font1 begin_bottom_legend_font2 begin_horizontal_scale begin_horizontal_scale_font begin_north_arrow begin_item begin_colorscale_indicator</Keywords>
           <Keywords name="Folder-">end_view <end_bottom_legend end_bottom_legend_font1 end_bottom_legend_font2 end_horizontal_scale end_horizontal_scale_font end_north_arrow end_item end_colorscale_indicator</Keywords>
           <Keywords name="Operators">. < ></Keywords>
           <Keywords name="Comment"></Keywords>
           <Keywords name="Words1"></Keywords>
           <Keywords name="Words2">WINDOW_HEIGHT WINDOW_WIDTH WINDOW_STATE VISIBLE X_ORIGIN Y_ORIGIN HEIGHT WIDTH SCALE_SIZE SCALE_TYPE DISPLAY_VIEW BACK_COLOR_R BACK_COLOR_G BACK_COLOR_B EXTENDS_WIDTH EXTENDS_HEIGTH EXTENDS_X EXTENDS_Y ARROW_TYPE BOTTOM_LEGEND_SHOW BOTTOM_LEGEND_SHOWLINE1 BOTTOM_LEGEND_LINE1 BOTTOM_LEGEND_SHOWLINE2 BOTTOM_LEGEND_LINE2 BOTTOM_LEGEND_COLOR_R BOTTOM_LEGEND_COLOR_G BOTTOM_LEGEND_COLOR_B BOTTOM_LEGEND_SHOWTIME TEXT_COLOR_R TEXT_COLOR_G TEXT_COLOR_B FILENAME ITEM_TYPE ITEM_CHECKED SHOW_SCALE HOR_ALIGNMENT VERT_ALIGNMENT TEXT_FONT TEXT_STYLE LEGEND_TEXT DIVISIONS SCALE_HORIZONTAL SCALE_ALPHA SCALE_BACKCOLOR_R SCALE_BACKCOLOR_G SCALE_BACKCOLOR_B POINT_SIZE POINT_SIZE_BY_PIXEL SHOW_POINT_NAMES OFF_SET_POINTS OFF_SET_X OFF_SET_Y MAX_VALUE MIN_VALUE TRANSPARENCY BELOW_LIMIT ABOVE_LIMIT BTN_CONSTANT_R BTN_CONSTANT_G BTN_CONSTANT_B BTN_BELOW_R BTN_BELOW_G BTN_BELOW_B BTN_ABOVE_R BTN_ABOVE_G BTN_ABOVE_B REPRESENTATION COLOR_SCALE_IDX LOGARITHMIC MIN_WIDTH MAX_WIDTH MIN_AREA X_STEP Y_STEP BTN_COLOR_R BTN_COLOR_G BTN_COLOR_B FILE SHOWCOLOR RENDER_COLOR INTERPOLATE_COLOR SHOWGRID BELOW_FILL</Keywords>
           <Keywords name="Words3"></Keywords>
           <Keywords name="Words4"></Keywords>
       </KeywordLists>
       <Styles>
           <WordsStyle name="DEFAULT" styleID="11" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="FOLDEROPEN" styleID="12" fgColor="008080" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="FOLDERCLOSE" styleID="13" fgColor="008080" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="KEYWORD1" styleID="5" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="1" />
           <WordsStyle name="KEYWORD2" styleID="6" fgColor="808080" bgColor="FFFFFF" fontName="" fontStyle="2" />
           <WordsStyle name="KEYWORD3" styleID="7" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="KEYWORD4" styleID="8" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT" styleID="1" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="COMMENT LINE" styleID="2" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="NUMBER" styleID="4" fgColor="FF8080" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="OPERATOR" styleID="10" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER1" styleID="14" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER2" styleID="15" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
           <WordsStyle name="DELIMINER3" styleID="16" fgColor="000000" bgColor="FFFFFF" fontName="" fontStyle="0" />
       </Styles>
   </UserLang>
 </NotepadPlus>
```
