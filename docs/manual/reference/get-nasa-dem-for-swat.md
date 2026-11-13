# Get NASA DEM for SWAT

This is a 3 step task and for simplification a structure as the one in \\\\Datacenter\\Albufeiras\_Eutrofizacao\\ComoFazer\\EstruturaDirectoriasAlbufeiras\\1.PreparaInputs\\3.SWAT\\1.MDT is advised

## Get the DEM from NASA

Obtain the DEM data from the [NASA site](http://edcsns17.cr.usgs.gov/EarthExplorer/)

1. On number 1 select "Digital Elevation" and check "SRTM"
2. On Number 2 select your coordinates (plus sign for two cordinate points)
3. Click Search
4. Wait for the Download to appear and save the results in DTED format.
5. Because of the changed web site the DTM is now no longer cutted in the area that you select, but comes in large areas. It is needed to cut within your area.
6. Open it in ArcGIS and save it as a raster,
7. Create a poligon in ArcCatalog and open it in ArcGIS
8. Change the poligon to the watershed area.
9. Use Spatial Analyst to reduce the extent of the DTM to the polygon.

Ask David For Help if you do not know to do these steps. A new video release should appear for the new steps because of the web site change.

## Retirar pixeis com NODATA no DEM (em PT - mini manual já estava escrito)

Por vezes o DEM pode vir sem dados em alguns pixeis. Este processo é para os retirar.

1. Colocar o DEM numa pasta em que o caminho seja simples e sem espaços no nome.
2. Ir a Programas -> ArcGIS -> ArcInfo Workstation -> ArcTools.
3. Aparece uma janela DOS e em cima, o menu ArcTools. Escolher >Grid Tools e carregar >ok.
4. Aparece mais uma janela DOS e a barra de menus GridTools.
5. No menu >ArcTools escolher >WorkSpace e definir a directoria onde está a GRID que se quer trabalhar.
6. No menu >Analysis escolher >Expression.
7. Aparece a janela GRID Expression Builder. A grid já deve aparecer onde diz >Input Grids.
8. Definir o nome da grid de output na caixa >Output Grid:.
9. Seleccionar NODATA conversion na caixa >Tool Categories.
10. Seleccionar NODATA to Value. Aparece uma janela para preencher- NODATA Conversion Funtion.
11. Na caixa >Input Grid: escrever o nome da grid de entrada. Na caixa >value: colocar o número zero (em número) e carregar >OK.
12. O comando foi escrito na caixa >Expression e basta carregar no botão >Apply para executar.
13. Se no primeiro passo foi copiada a grid original para uma pasta mais simples, copiar a grid gerada deste processo para a pasta de trabalho (não esquecer da pasta info)

## Change the coordinate system from Geographic to Militares Portuguesas

Check second part of the video ComoObterMDTParaSWAT.swf at \\\\\\Datacenter\\Albufeiras\_Eutrofizacao\\Produtos\\Apresentacoes\\Videos
