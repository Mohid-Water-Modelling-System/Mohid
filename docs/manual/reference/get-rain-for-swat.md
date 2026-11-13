# Get Rain For SWAT

Check video ComoObterPrecipitcaoParaSWAT.swf in \\\\datacenter\\Albufeiras\_Eutrofizacao\\Produtos\\Apresentacoes\\Videos

1. Use MOHID GIS or ArcGIS with the rivers shape to check all the rain station in your watershed.
2. Go to snirh.pt and collect the daily rain for all of them. If you have only one or two stations inside the watershed, get stations outside but near, because you may have periods with no data and need to fill.
3. Use macro excel file OrdenaPrecip.xls in \\\\Datacenter\\Albufeiras\_Eutrofizacao\\ComoFazer\\1.PreparaInputs\\1.Dados\\Meteorologia to get all the precipitation files side by side and with no day gaps.
4. Fill the gaps in the stations: i) if the gap is only one day or two use the average of the last existing days; ii) if gaps are bigger use the values from a station nerby and that is on the same "Isoieta" (same precipitation lines) - see Portuguese Isoietas (PRET) shape in \\\\datacenter\\Albufeiras\_Eutrofizacao\\Produtos\\ProjectosGIS\\ShapesBase\\AtlasAmbiente\_Isoietas
5. Save each station in a separate file (.dbf) with the format as in \\\\datacenter\\Albufeiras\_Eutrofizacao\\ComoFazer\\1.PreparaInputs\\3.SWAT\\ExemploAplicacao\\data\\ValadaPCP.dbf.
