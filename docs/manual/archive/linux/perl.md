# Perl

Perl is a powerful and versatile scripting language that comes in handy when bot-working is needed on any operating system.

## Installing

To use the perl scripting language, first you need to install it on your operating system.

- To install packages type:

```
>perl -MCPAN -e shell
CPAN> install packagename
CPAN> force install packagename
```

(ex: install Date::Calc). **Warning**: under windows use a .NET environment command prompt!!

- To run a perl script type:

```
>perl scriptfile
```

(ex: perl test.pl).

### Installing for windows

If you are one of the unlucky Visual Studio 2005 owners, then you'd better install the [CamelPack](http://sourceforge.net/project/showfiles.php?group_id=158775&package_id=178025) or else you won't be able to install perl modules from CPAN. Microsoft did there a hard blow against the Perl community.

### Common extensions

Here are the common extensions used:

```
use Net::FTP;
use Net::SMTP;
use IO::File;
use Date::Calc;
use NetCDF;
```

### PPM

The Perl Package Manager is here to help installing and managing Perl packages. The greatest asset of Perls longevity is its package repository with over 6500 packages! Here are some of the commands that best explain how to use ppm:

```
> ppm list
```

will list installed packages.

```
> ppm help
> ppm help command
```

will return help manual.

```
> ppm search pattern
> ppm search *
```

will return a numbered list of downloadable packages matching *pattern*.

```
> ppm describe pkg
> ppm describe module
> ppm describe url
> ppm describe number
```

will describe the selected package. Note that *number* is the package number taken from the last search.

```
> ppm install pkg
> ppm install module
> ppm install url
> ppm install number
```

will install the desired package. Note that *number* is the package number taken from the last search.

## Variables, arrays and hashes

### Variables

### Arrays

### Hashes

Here's a sample hash construction

```
   %hash = (
       key1 => 'value1',
       key2 => 'value2',
       key3 => 'value3',
   );
```

Here's a simple hash entry definition

```
$hash{$key} = $value;
```

## Flow control

### For

### While

```
while ( ($key, $value) = each(%$nc_ref) ) {
   print "$key => $value \n";
}
```

### Foreach

```
foreach $thing (@list){
    $machine = $thing;
    $filesize = $thingy;
    $filename = $thingabob;
    write;
}
```

## [Regexp](regexp.md)

### matching

Returns true if *pattern* is found in *$sample*

```
$sample =~ m/pattern/
```

### Substituting

Substitutes the *pattern* with *newtext* in every occurence of string *sample* and writes the changes into *sample*.

```
$sample =~ s/pattern/newtext/g
```

### Translate

Translates an ascii letter for another. It applies to patterns as well.

```
$sample =~ tr/,\./; /g
```

The above examples changes the commas(,) for semi-colons(;) and points(.) for spaces( ).

## [One-Liners](perl-one-liners.md)

Perl is mighty powerful from the shell command line, simply by invoking [one-liners](perl-one-liners.md) and pipes.

## Formatting text output

This snippet typifies how to format text in perl.

- **\>** justifies right
- **<** justifies left
- **|** justifies center

```
format =
@>>>>>>>>>>>>>> @>>>>>>>>>> @>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
$machine, $filesize, $filename
.
```

```
foreach $thing (@list){
    $machine = $$thing[0];
    $filesize = $$thing[1];
    $filename = $$thing[2];
    write;
}
```

## Sample files

### Rename files with a regexp

```
#!/usr/local/bin/perl
#
# Usage: rename perlexpr [files]

($regexp = shift @ARGV) || die "Usage:  rename perlexpr [filenames]\n";

if (!@ARGV) {
   @ARGV = <STDIN>;
   chomp(@ARGV);
}

foreach $_ (@ARGV) {
   $old_name = $_;
   eval $regexp;
   die $@ if $@;
   rename($old_name, $_) unless $old_name eq $_;
}

exit(0);
```

Example:

```
> perl perlrename.pl "s/_/-/g" test_file_.txt
> dir
>      test-file-.txt
```

Now install it in c:/Perl/bin and it'll be universally accessible

```
> perlrename.pl "s/_/-/g" test_file_.txt
> for /F %i in ('dir /B *.eps') do @perlrename.pl "s/_/-/g" %i
```

### Example of a script skeleton

`[perl,Y]`

```
#!/usr/bin/perl -w
use strict;

##
## Usage Information
##

my $usage = <<EOF;

$0 <argument1> <argument2>
This script does this and needs <arguments>

EOF

##
## Parse Arguments
##

if ($#ARGV < 1) {
    print $usage;
    exit 1;
};

my $arg1 = shfit(@ARGV);
my $arg2 = shfit(@ARGV);

##
## Begin
##

$stat = system("ls *");
&DieNice ($stat, "ERROR on ls *");

##
## Subroutines
##

sub DieNice {

   my ($stat, $msg) = @_;
   die "$msg \n" if ($stat !=0);

};
```

### Manipulate [netcdf](netcdf.md) files

`[perl,Y]`

```
#!/usr/bin/perl -w
#
#netcdf rivers {
#dimensions:
#        river = 3 ;
#        time = UNLIMITED ; // (366 currently)
#variables:
#        int river(river) ;
#                river:title = "Number of river" ;
#                river:long_name = "1:Loire 2:Gironde 3:Adour" ;
#        float time(time) ;
#                time:units = "seconds since 2004-01-01 00:00:00" ;
#                time:title = "Time" ;
#        float runoff(time, river) ;
#                runoff:long_name = "Runoff of the river" ;
#                runoff:units = "m3s-1" ;
#        float temperature(time, river) ;
#                temperature:long_name = "Temperature" ;
#                temperature:units = ".C" ;
#        float salinity(time, river) ;
#                salinity:long_name = "Salinity" ;
#                salinity:units = "PSU" ;
#
#// global attributes:
#                :title = "Runoff data" ;
#                :institution = "GIP MERCATOR OCEAN" ;
#                :references = "http://www.mercator-ocean.fr" ;
#}
#---------------------------------

use NetCDF;
use strict;

# Input netcdf filename
my $infile = "rivers.nc";

#Two dimensions: time(366) and rivers(3). Stride one river at a time.
my %riv = (
    Loire => 0,
    Gironde => 1,
    Adour => 2,
);

# Choose between Loire, Gironde or Adour
my $key = 'Gironde';
print "Doing river $key\n\n";

my $len=366;
my @start = (0, 0);
my @count = ($len, 3);
my @slice =();
for (my $i=0; $i<$len; $i++) {
    push (@slice, $i*3 + $riv{$key});
}
print "@slice\n";

# Open netcdf file in read-only mode
my $ncid = NetCDF::open($infile, NetCDF::NOWRITE);
print "ncid: $ncid\n";

# Dimensions
my $d_rivid = NetCDF::dimid($ncid, 'river');
print "d_rivid: $d_rivid\n";

my $d_timeid = NetCDF::dimid($ncid, 'time');
print "d_timeid: $d_timeid\n";

# Variables
my @river = ();
my $rivid = readvar($ncid, 'river', (0), (3), \@river);
print "@river\n\n";

my @time = ();
my $timeid = readvar($ncid, 'time', (0), (366), \@time);
#print "@time\n\n";

my @runoff = ();
my $runoffid = readvar($ncid, 'runoff', \@start, \@count, \@runoff);
my @runoff_s = @runoff[@slice];
print "@runoff_s\n\n";

my @temperature = ();
my $temperatureid = readvar($ncid, 'temperature', \@start, \@count, \@temperature);
my @temperature_s = @temperature[@slice];
print "@temperature_s\n\n";

my @salinity = ();
my $salinityid = readvar($ncid, 'salinity', \@start, \@count, \@salinity);
my @salinity_s = @salinity[@slice];
print "@salinity_s\n\n";

# Creating time in days variables
my @timeday = ();
for (my $i=0; $i<$len; $i++) {
    push (@timeday,$time[$i]/86400);
}
print "@timeday\n\n";

# Close netcdf file
NetCDF::close($ncid);

# Write the TimeSeries dat file
writetimeserie (

    $key,
    \@timeday,
    \@runoff_s,
    \@temperature_s,
    \@salinity_s,

);

sub readvar {
# syntax:
#       my var = "";
#       my varid = readvar($ncid, $varname, \@start, \@count, \@var);

    my $ncid = shift;
    my $varname = shift;
    my $start = shift;
    my $count = shift;
    my $var_ref = shift;

    my $varid = NetCDF::varid($ncid, $varname);
    print "${varname}id: $varid\n";
    my @tempvar = (0,0);
    my $status = NetCDF::varget($ncid, $varid, $start, $count, \@tempvar);
    print "Status: $status\n";
    #print "@tempvar\n";
    #Copy array
    @$var_ref = @tempvar;

    return $varid;

}

sub writetimeserie {
# Syntax:
#   writetimeserie($rivername, \@timedays, \@runoff,  \@temperature, \@salinity);

    my $rivername = shift;
    my $timev = shift;
    my $runv = shift;
    my $tempv = shift;
    my $saltv = shift;

    print "$#$timev\n";

    #Open file to write
    $rivername.=".srw";
    open(FILE, ">$rivername");

    print FILE "TIME_UNITS  :   DAYS\n";
    print FILE "SERIE_INITIAL_DATA  :   2004    1   1   0   0   0\n\n";

    my $time = 'Tempo';
    my $flow = ' Flow';
    my $temp = '    T';
    my $salt = '    S';
    print FILE "$time $flow $temp $salt\n";

    $time = ' dias';
    $flow = ' m3/s';
    $temp = '    C';
    $salt = '  psu';
    print FILE "$time $flow $temp $salt\n";

    print FILE "\n";

    $time = '    1';
    $flow = '    2';
    $temp = '    3';
    $salt = '    4';
    print FILE "$time $flow $temp $salt\n";

    print FILE "\n";

    $time = 'Tempo';
    $flow = ' Flow';
    $temp = '    T';
    $salt = '    S';
    print FILE "$time $flow $temp $salt\n";

    print FILE "<BeginTimeSerie>\n";

    $time = '99';
    $flow = 99;
    $temp = 99;
    $salt = 99;

format FILE =
@<<<< @###.##### @###.# @###.#
$time, $flow, $temp, $salt
.

    for (my $i=0; $i<$#$timev+1; $i++){

        $time = $$timev[$i];
        $flow = $$runv[$i];
        $temp = $$tempv[$i];
        $salt = $$saltv[$i];
        write FILE;

    }

    print FILE "<EndTimeSerie>\n";

    # Close file
    close(FILE);

}
```

### Substitute pattern in file

```
#!/usr/local/bin/perl
#This program substitutes the pattern "dd-mm-yyyy hh:mm:ss" from
#argument file with "dd, mm, yyyy, hh, mm, ss".
#The new text is then written in an output file.

open FILE, "< @ARGV[0]";
@text=<FILE>;
close FILE;

foreach $line (@text)
{
      $line =~ s/(\S{1,2})-(\S{1,2})-(\S{4}) (\S{2}):(\S{2}):(\S{2})/$1, $2, $3, $4, $5, $6/g;
      push @newdata, $line;
}

open FILE, "> new_$ARGV[0]";
print FILE @newdata;
close FILE;
```

### Send an e-mail

Here's a script that allows one to send e-mails. `[perl,Y]`

```
use Net::SMTP;

$HOST = 'mail.ist.utl.pt';
$FROM = 'my_address';
$TO = 'his_address';
$SUBJECT = 'SMTP perl test';
@MESSAGE = (
               "JAPH\n",
		"Maretec team\n",
	    );
$smtp = Net::SMTP -> new(	Host => $HOST,
				Hello => $FROM,
				Timeout => 60,
				Debug => 1
			);
$smtp -> mail($FROM);
$smtp -> to($TO);
@mail_content = (	"Subject: ".$SUBJECT."\n",
			"From: $FROM\n",
			"To: $TO\n",
			"\n"
		);
push(@mail_content,@MESSAGE);
$smtp -> data(@mail_content);
$smtp -> quit();
```

### FTP access

Here's a script to download from a remote ftp server `[perl,Y]`

```
use Net::FTP;

$SERVER = 'your.server';
$USER = 'you';
$PASSW = 'yourpass';
$FTP_DIR = 'yourdir';

# Login and cd ( select binary mode ).
$ftp = Net::FTP->new("$SERVER") || die "Can't communicate with $SERVER!";
$ftp->login($USER, $PASSW) || die "Can't login $USER,$PASSW!";
$ftp->binary || die "Can't change mode!";
$ftp->cwd($FTP_DIR) || die "can't open $FTP_DIR";

# Get the files list of the ftp.
my @files = $ftp->ls();

# Download the files list.
foreach $file (@files){
	$ftp->get($file) || die "Can't download $file: ".$ftp->message."\n";
	print "Downloaded $file\n";
}
$ftp->quit();
```

### rename files

Here's a sample file that renames all the files containing *.f90* to *.F90* `[perl,Y]`

```
#!/usr/local/bin/perl
#
# Program to substitute $ren to $to for all files in the directory
#

use strict;

$ren = ".f90";
$to = ".F90";

opendir(DIR,".");
@text = readdir(DIR);
foreach $line (@text)
{
	if ($line =~ /$ren/)
	{
		@words = split(/ /, $line);
		foreach $word (@words)
		{
			if ($word =~ /$ren/)
			{
				$word =~ s/\n//;
				push (@files, $word);
			}
		}
	}
}
closedir(DIR);

foreach $file (@files)
{
	$file2 = $file;
	$file2 =~  s/$ren/$to/;
	rename($file, $file2) || die "Can't rename ".$file;
	print "Changed ".$file." to ".$file2."\n";
}
```

### ConvertToHDF5

Here's another sample file that I use to automatically [ConvertToHDF5](../../tools/converttohdf5.md) all the compressed [netcdf](netcdf.md) files in a given repository: `[perl,Y]`

```
#!/usr/local/bin/perl
#
#

use strict;

$DATE = '20060628';
$IST = 'ist_meteog-mercatorPsy2v2r1v_R';
$FOLDER = $IST . $DATE;
$GZIP = gzip;
$NC = ".nc";
$CONVERT = "ConvertToHDF5";

# 1 Uncompress them ###############

chdir($FOLDER) or die "Can't chdir to ".$DIR;
print 'Uncrunching ...\n';
@args = ($GZIP, "-d", "\*\.gz");
system(@args) == 0 or print "system @args failed: $?\n";

opendir(DIR, ".");
@text = readdir(DIR);
foreach $line (@text)
{
	if ($line =~ /$NC/)
	{
		@words = split(/ /, $line);
		foreach $word (@words)
		{
 			if ($word =~ /$ren/)
			{
				$word =~ s/\n//;
				push (@files, $word);
			}
		}
	}
}
print "Fetched ".@files." files.\n";
chdir("..");

# 2 Generate ConvertToHDF5Action ###############

$action = 'ConvertToHDF5Action.dat';
open(FILE, ">$action");
my $TXT = << EndOfTXT;

<begin_file>
ACTION                    : CONVERT MERCATOR FORMAT
OUTPUTFILENAME            : $FOLDER/MercatorR$DATE.hdf5
OUTPUT_GRID_FILENAME      : $FOLDER/MercatorGridR$DATE.dat
INPUT_GRID_FILENAME       : GridFiles/ist_meteog-gridT.nc
INPUT_GRID_FILENAME_U     : GridFiles/ist_meteog-gridU.nc
INPUT_GRID_FILENAME_V     : GridFiles/ist_meteog-gridV.nc
OUTPUT_GEOMETRY_FILENAME  : $FOLDER/MercatorGeometryR$DATE.dat
<<begin_input_files>>

EndOfTXT;

print FILE $TXT;

foreach $file (@files)
{
	print FILE $FOLDER."\/".$file."\n";
}
print FILE "<<end_input_files>>\n";
print FILE "<end_file>\n";
close(FILE);

# 3 Run ConvertToHDF5 ###############

print "Converting to HDF5...\n";
@args = ($CONVERT);
system(@args) == 0 or die "system @args failed: $?";
print "Done ".@files." files.\n";

# 4 Compress them back again ####################

chdir($FOLDER) or die "Can't chdir to ".$DIR;
print "Crunching back up...\n";
@args = ($GZIP, "\*\.nc");
system(@args) == 0 or die "system @args failed: $?";
print "Finished ".@files." files.\n";
```

## References

1. The wikipedia [article](http://en.wikipedia.org/wiki/Perl_programming_language).
2. A popular [tutorial](http://www.comp.leeds.ac.uk/Perl/start.html).
