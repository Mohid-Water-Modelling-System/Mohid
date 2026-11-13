# Kml

[KML](kml.md) is the Keyhole Markup Language. It's an xml-based language that allows to create content for Google Earth and Google Maps.

It allows to create viewpoints, placemarks, paths, shapes, image overlays and 3D buildings. It also allows to embed customized [html](../web/html.md) code in the placemarks.

## Sample kml file

### Placemarks with custom html code

[Here's](http://maps.google.com/maps?f=q&q=http%3A%2F%2Fnfist.ist.utl.pt%2F%7Eguillaume%2FMohidBeachcams_NetL.kml&btnG=Search+Maps) a simple **kml** file that adds a placemark, a viewpoint for it, and embeds some custom html code:

```
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.2">

<Document>

 <name>MohidBeachCams.kml</name>

 <visibility>1</visibility>

 <Folder>

   <name>Portugal Beach Cams</name>

   <visibility>1</visibility>

   <description>
    <![CDATA[
     Live data from <a href="http://cams.beachcam.pt/">Beachcam_PT</a>
     and <a href="http://www.surftotal.com/webcams2/Default.htm">SurfTotal webcams</a>.
     Best viewed when using GE built-in web browser.
     Tools -> Options -> Preferences -> Navigation bar -> Uncheck "Show web results in external  browser".
   ]]>
   </description>

   <LookAt>
     <longitude>-9.285977261935919</longitude>
     <latitude>38.7380813109443</latitude>
     <range>47801.2002729008</range>
     <tilt>2.43637162447063e-011</tilt>
     <heading>0.01425072853298695</heading>
   </LookAt>

   <Placemark>

     <name>Praia do Guincho</name>

     <LookAt>
       <longitude>-9.470984593789977</longitude>
       <latitude>38.73175408910536</latitude>
       <range>1421.36386586457</range>
       <tilt>2.897405864852683e-010</tilt>
       <heading>-0.08870325629105799</heading>
     </LookAt>

     <styleUrl>#camneutra</styleUrl>

     <Point>
       <coordinates>-9.472970249051258,38.73236056451183,0</coordinates>
     </Point>

     <description>
        <![CDATA[
         Looking from Praia do Guincho.
         <img src="http://cams.beachcam.pt//camaxis/guincho_00001.jpg">
         Surf report <a href="http://cams.beachcam.pt/praias_beachcams.php?id=18">here</a>.
         Live cam <a href="http://www.surftotal.com/webcams2/guincho_cam.asp?id=25">here</a>.
        ]]>
     </description>

   </Placemark>

 </Folder>

</Document>

</kml>
```

### Network links

Network links can be very useful when providing early releases of beta services. This way, the developper may continue to upgrade the production version of the **kml** file, by providing to the users a file that contains a link to it.

[Here's](http://maps.google.com/maps?f=q&q=http%3A%2F%2Fnfist.ist.utl.pt%2F%7Eguillaume%2FMohidBeachcams_NetL.kml&btnG=Search+Maps) a simple kml file that adds a network link to another **kml** file:

```
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.2">
<NetworkLink>
       <name>MohidBeachcams</name>
       <open>1</open>
       <flyToView>1</flyToView>
       <Url>
               <href>http://nfist.ist.utl.pt/~guillaume/BeachCams_PT.kml</href>
               <refreshMode>onInterval</refreshMode>
               <refreshInterval>50</refreshInterval>
               <viewRefreshMode>onRegion</viewRefreshMode>
       </Url>
</NetworkLink>
</kml>
```

### Image overlays

### Regions

### Time

## Kml2Mohid

[Kml2Mohid](../../development/programming/kml2mohid.md) is a little [ruby](ruby.md) program that converts [Google Earth](../web/google-earth.md) **kml** files into Mohid GIS data files.

## External References

- [Kml reference documentation](http://code.google.com/apis/kml/documentation/index.html)
- [MohidBeachCams](http://maps.google.com/maps?f=q&q=http%3A%2F%2Fnfist.ist.utl.pt%2F%7Eguillaume%2FMohidBeachcams_NetL.kml&btnG=Search+Maps)
- [Google Earth](http://earth.google.com)
- [Google Maps](http://maps.google.com)
- [Google Mapplets](http://maps.google.com/maps/mm?mapprev=1)
