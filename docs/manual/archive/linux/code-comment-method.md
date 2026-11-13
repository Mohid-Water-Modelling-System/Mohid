# Code Comment Method

A good practice for any language programming is to comment abundantly the code. This way, things will be easier for the coder, or any code reader, to decode in the future the meaning of the code. However, each coder has its own standards to comment code; indeed as some coders are very succint, others are excessively prolific, while others remain cryptic at best.

Hence, a standard for commenting code would prove useful. To do so, one needs to establish rules for code commenting, and Mohid does have a set of rules. However, these are quite flexible, and most of all, the implemented comments in the code is ill-formed nearly everywhere.

The advantages to establish a strict set of rules are immense! The greatest one of them would be to automatically generate code documentation!

Hence, basing ourselves in [XML](http://en.wikipedia.org/wiki/Xml)/[XSD](http://en.wikipedia.org/wiki/Xsd)/[XSLT](http://en.wikipedia.org/wiki/Xsl) technology we will outline the principles for a sound code commenting practice. The method is currently applied to comment perl scripts, thus only the respective xml and schema are presented here.

## Commenting code using the XML method

To comment code one needs to define for each program:

- The program name
- A header that contains:
    -   A description of the program
    -   The program user's syntax
    -   The program inputs
    -   The program outputs
- The code that contains:
    -   The generic comments
    -   The subroutines, which in turn contain:
        -   A description
        -   The inputs
        -   The outputs

All of these directives are well-accepted and of common sense. However we decide to be very strict to them and decide to use XML to write the comments in the code.

### Sample perl code

Here is a sample perl script that uses the commenting XML method:

```
#!/usr/local/bin/perl
#<PROGRAM_DOCUMENT>
#
#<TYPE>Perl</TYPE>
#<PROJECT>Operational Model Automated Scripts</PROJECT>
#<NAME>ExtractComments</NAME>
#
#<HEADER>
#
# <DESCRIPTION>
# 	Program that extracts the comments from the perl scripts and outputs in
# 	an XML format defined by a XSD schema.
# 	!As a rule the first line is ignored!
# </DESCRIPTION>
#
# <SYNTAX>
#	$perl ExtractComments.pl < perlscript.pl > _perlscript.xml\n
#	$more perlscript.pl | perl ExtractComments.pl > _perlscript.xml.
# </SYNTAX>
#
# <INPUT>
#	<STDIN>void</STDIN>
#      <FILES>
#          <FILE>perlscript.pl</FILE>
#      </FILES>
# </INPUT>
#
# <OUTPUT>
#	<STDOUT>
#           Comments in code.
#	</STDOUT>
# </OUTPUT>
#
#</HEADER>
#
#<_CODE>
use IO::File;

if (!@ARGV)
{
   @ARGV = <STDIN>;
   chomp(@ARGV);
}
#What if this happens?
while (!($ARGV[0] =~ /<PROGRAM_DOCUMENT>/i)) #Are we at the beginning of the document yet?
{
	shift @ARGV; #No? Then keep looking.
}
foreach $line (@ARGV)
{
	$_= substr($line,0,1); #Isolate the first character
	if ($_ =~ /#/) #Is it a comment char?
	{
 		$newline = substr($line,1); #Extract the line without the "#" character.
	 	push (@outtext, "$newline\n"); #and push it to the output.
	}
}

print "@outtext\n"; #output to <STDOUT>
# <SUBROUTINE>
# <NAME>openlog</NAME>
# <DESCRIPTION>
# 	This routine creates or opens a log file and appends a timestamp.
# </DESCRIPTION>
# <INPUT>
# 	logfilename.
# </INPUT>
# <OUTPUT>
# 	$LOG IO::File handler.
# </OUTPUT>
# </SUBROUTINE>
sub openlog
{
	$timestamp=localtime();
	$LOG = new IO::File ">>$_[0]";
	print $LOG "\n############# $timestamp $_[0] ######################\n\n";
	return $LOG;
}
#</_CODE>
#</PROGRAM_DOCUMENT>
```

## Extracting the comments from the code

At first, it may seem a bit absurd to use XML to comment code as it turns it a lot heavier than before. However, the perl sample script extracts the comments from the code. Here's the command:

```
$perl ExtractComments.pl < ExtractComments.pl
```

### Comments from the sample perl code

And here's the output:

```
<PROGRAM_DOCUMENT>

<TYPE>Perl</TYPE>
<PROJECT>Operational Model Automated Scripts</PROJECT>
<NAME>ExtractComments</NAME>

<HEADER>

 <DESCRIPTION>
       Program that extracts the comments from the perl scripts and outputs in
       an XML format defined by a XSD schema.
       !As a rule the first line is ignored!
 </DESCRIPTION>

 <SYNTAX>
       $perl ExtractComments.pl < perlscript.pl > _perlscript.xml\n
       $more perlscript.pl | perl ExtractComments.pl > _perlscript.xml.
 </SYNTAX>

 <INPUT>
	<STDIN>void</STDIN>
      <FILES>
          <FILE>perlscript.pl</FILE>
      </FILES>
 </INPUT>

 <OUTPUT>
	<STDOUT>
           Comments in code.
	</STDOUT>
 </OUTPUT>

</HEADER>

<_CODE>
What if this happens?
 <SUBROUTINE>
 <NAME>openlog</NAME>
 <DESCRIPTION>
       This routine creates or opens a log file and appends a timestamp.
 </DESCRIPTION>
 <INPUT>
       logfilename.
 </INPUT>
 <OUTPUT>
       $LOG IO::File handler.
 </OUTPUT>
 </SUBROUTINE>
</_CODE>
</PROGRAM_DOCUMENT>
```

Neat, huh?

#### The exception rule

The only exception rule is that the `<, >` and `&` characters cannot be displayed as is in XML format. The recommended alternative way to represent them is:

```
<ampersand>lt; , <ampersand>gt; and <ampersand>amp;
```

## Validating the comments from the code

Indeed it's very interesting as we get "ready to use" XML code to produce a program documentation! However, well-formed xml is required to produce documentation. Furthermore, an XML Schema is also required if one wants to use a XSLT stylesheet to port the XML to HTML or PDF. Thus, validate the code comments is required in order to guarantee that code comments are well-formed and conform with the rules.

### The Schema

Here is the XML Schema (XSD) that conforms to the example given above:

```
<?xml version="1.0"?>

<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">

	<-- Root Level -->
	<xsd:element name="PROGRAM_DOCUMENT" type="PROGRAM_DOCUMENTtype"/>

	<-- Zeroth Level -->
	<-- 1 type to be defined here: PROGRAM_DOCUMENT -->
	<xsd:complexType name="PROGRAM_DOCUMENTtype">
		<xsd:sequence>
			<xsd:element name="TYPE" type="TYPEtype"/>
			<xsd:element name="PROJECT" type="xsd:string" minOccurs="0"/>
			<xsd:element name="NAME" type="xsd:string"/>
			<xsd:element name="HEADER" type="HEADERtype"/>
			<xsd:element name="CODE" type="CODEtype"/>
		</xsd:sequence>
	</xsd:complexType>

	<-- First Level -->
	<-- 3 types to be defined here: TYPE, HEADER, CODE -->
	<xsd:simpleType name="TYPEtype">
		<xsd:restriction base="xsd:string">
			<xsd:enumeration value="Perl"/>
			<xsd:enumeration value="Fortran77"/>
			<xsd:enumeration value="Fortran90"/>
			<xsd:enumeration value="Fortran95"/>
	 		<xsd:enumeration value="VB"/>
			<xsd:enumeration value="Php"/>
			<xsd:enumeration value="Cgi"/>
			<xsd:enumeration value="Python"/>
		</xsd:restriction>
	</xsd:simpleType>

	<xsd:complexType name="HEADERtype">
		<xsd:sequence>
			<xsd:element name="DESCRIPTION" type="xsd:string"/>
			<xsd:element name="SYNTAX" type="xsd:string"/>
			<xsd:element name="INPUT" type="INPUTtype"/>
			<xsd:element name="OUTPUT" type="OUTPUTtype"/>
		</xsd:sequence>
	</xsd:complexType>

	<xsd:complexType name="CODEtype" mixed="true">
 		<xsd:sequence>
			<xsd:element name="SUBROUTINE" type="SUBROUTINEtype" minOccurs="0"  maxOccurs="unbounded"/>
		</xsd:sequence>
	</xsd:complexType>

	<-- Second Level -->
	<-- 3 types to be defined here: INPUT, OUTPUT, SUBROUTINE -->
	<xsd:complexType name="INPUTtype" mixed="true">
		<xsd:sequence>
			<xsd:element name="STDIN" type="STDtype" minOccurs="0"/>
			<xsd:element name="FILES" type="FILEStype" minOccurs="0"/>
		</xsd:sequence>
	</xsd:complexType>

	<xsd:complexType name="OUTPUTtype" mixed="true">
		<xsd:sequence>
			<xsd:element name="STDOUT" type="STDtype" minOccurs="0"/>
			<xsd:element name="FILES" type="FILEStype" minOccurs="0"/>
		</xsd:sequence>
	</xsd:complexType>

	<xsd:complexType name="SUBROUTINEtype" mixed="true">
		<xsd:sequence>
			<xsd:element name="NAME" type="xsd:string"/>
			<xsd:element name="DESCRIPTION" type="xsd:string"/>
			<xsd:element name="INPUT" type="xsd:string"/>
			<xsd:element name="OUTPUT" type="xsd:string"/>
	 	</xsd:sequence>
	</xsd:complexType>

	<-- Third Level -->
	<-- 2 types to be defined here: STDtype, FILEStype -->
	<xsd:complexType name="STDtype" mixed="true">
		<xsd:sequence>
			<xsd:element name="EXAMPLE" type="xsd:string" minOccurs="0"/>
		</xsd:sequence>
	</xsd:complexType>

	<xsd:complexType name="FILEStype">
		<xsd:sequence>
			<xsd:element name="FILE" type="xsd:string" minOccurs="0"/>
		</xsd:sequence>
	</xsd:complexType>

</xsd:schema>
```

### The validate action

Now that we have the comments in XML format and an XSD file to validate them against, we'll be requiring a command line utility that performs the validation. Here we will use the free java program [validateXML.jar](http://www.mohid.com/hydrogroup/data_sources/validateXML.jar) and some shell scripting (piping) to perform the trick:

#### The XML header

First we need to give the XML comments a header that links the XML to the XSD. We use a little [perl one-liner](perl.md) for that:

```
$perl -p -e "s/<PROGRAM_DOCUMENT>/<?xml version = \"1.0\"?>\n\n<PROGRAM_DOCUMENT xmlns:xsi =
\"http:\/\/www.w3.org\/2001\/XMLSchema-instance\" xsi:noNamespaceSchemaLocation =
\"http:\/\/www.mohid.com\/hydrogroup\/data_sources\/metacode.xsd\">/i;" < ExtractComments.xml
```

Here's an extract of the output:

```
<?xml version = "1.0"?>

<PROGRAM_DOCUMENT xmlns:xsi = "http://www.w3.org/2001/XMLSchema-instance"
xsi:noNamespaceSchemaLocation = "http://www.mohid.com/hydrogroup/data_sources/metacode.xsd">

<TYPE>Perl</TYPE>
<PROJECT>Operational Model Automated Scripts</PROJECT>
<NAME>ExtractComments</NAME>
...
```

#### The validation check

Next we can check to see if the XML is well-formed and conform to its schema with [validateXML.jar](http://www.mohid.com/hydrogroup/data_sources/validateXML.jar):

```
$java -jar validateXML.jar ExtractComments_new.xml
```

Here's the output if the XML is valid:

```
Validating "ExtractComments_new.xml":
root element = "":PROGRAM_DOCUMENT
OK
```

## Producing a document

If this method is rigorously applied and no program may be checked in the local documents repository until it passes the schema validation test, then structured coherent comments are readily available to be ported in a more reader-friendly format. The best way is to produce a XSLT stylesheet that transforms the XML into, say, a pdf format or openoffice, or word. This has [yet](http://www.mohid.com/hydrogroup/codecomments/ExtractComments.pl.xml) to be implemented in this perl example. Another very interesting idea would be to post the code comments to a blog or a wiki or a web-page(via mail).

## Implementation

Currently, the perl commenting method files can be found in SourceSafe in the OtherTools database under the folder Perl/Meta. Also, an experimental [web-based code repository](http://www.mohid.com/hydrogroup/codecomments.xml) is being tested.

## External References

- [W3C XML.](http://www.w3.org/XML/)
- [W3C XML Schema.](http://www.w3.org/XML/Schema)
- [W3School XSD](http://www.w3schools.com/schema/default.asp) is a popular tutorial for XSD.
- [XmlStartlet](http://xmlstar.sourceforge.net/) is a very useful suite of command line tools.
- [validateXML.jar](http://www.mohid.com/hydrogroup/data_sources/validateXML.jar) is the java program used to validate the XML against the XSD.
