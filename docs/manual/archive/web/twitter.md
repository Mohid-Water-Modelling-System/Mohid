# Twitter

[Twitter](http://www.twitter.com) is a social short messaging system from the web to the outside world. It can be used to send messages to friends. A practical utility is to have a bot-twitter-user that broadcasts messages from scheduled scripts. It could confirm the status of an ongoing process or job with an sms.

## Using twitter with [curl](../linux/curl.md)

```
curl -u mohid1:***** -d status="Having fun with curl" http://twitter.com/statuses/update.xml
```

## Mohid Twitter project

Make sure you have accessible from your system path the following GNU utilities:

1. tee
2. grep
3. tail
4. curl

Save this code as a batchfile:

```
REM MohidTwitter.bat
REM Usage syntax:
REM > MohidWater.exe | tee Run1.log
REM > Call MohidTwitter.bat Run1.log
@echo off

REM Please edit the username and password value
set username=MyUsername
set password=MyPass

type %1 | grep MOHID | grep "[^-]$" | tail -n 1 > _tmp.txt
set /P mess=<_tmp.txt
echo %mess%
curl -u %username%:%password% -d status="%mess%" "http://twitter.com/statuses/update.xml"
del _tmp.txt

@echo on
```

## Bot agents

- [mohid1](http://www.twitter.com/mohid1): Pre-op model;
- [weatherlisbon](http://www.twitter.com/weatherlisbon): BBC lisbon weather forecast;
- [marvin\_bot](http://www.twitter.com/marvin_bot): Marvin, the paranoid android;
- [euromillions](http://www.twitter.com/euromillions): Every friday night you get the euromillions results;
- [weatherhannover](http://www.twitter.com/weatherhannover): BBC hannover weather forecast;
- [Biosfera](http://www.twitter.com/biosfera): Meteo-IST;

## Maretec agents

- [Guillaume R.](http://www.twitter.com/guillaume7)
- [Bartolomeu B.](http://www.twitter.com/bartb_pt)

| Channel Title: |  |
| --- | --- |
|  |

## Twitter bot script

### bash example

Here's a sample bash script that uses curl to extract a weather forecast from an rss feed and twit it:

`[bash,Y]`

```
#! /bin/sh

#This line is important if you are to put this script as a cronjob
#Goto the folder where the bash script exists.
here=/home/[user]/scripts
cd $here

#BBC Lisbon weather id
id=0048

#BBC weather RSS feed address
feed="http://feeds.bbc.co.uk/weather/feeds/rss/5day/world/${id}.xml"

#City
city=lisbon

#temporary file
file="weather${city}.txt"

#Weather twitter bot
twitbot=weatherlisbon:w*****r

#Read the RSS feed and filter it
curl $feed 2>&1 | grep 'title' | tail -n 1 | perl -wlne'm/title>(.*)<\/title/i && print $1' |
sed -e  "s/\&\#xB0;//g" > $file
# use this line instead, to filter out the Fahrenheit degrees, as well as the leading "the forecast for"
#sed -e "s/\&\#xB0;//g" | sed -e "s/ [\(][0-9]*F[\)]//g" | sed -e "s/space:\+/ /g" | cut -d ' ' -f4-21 > $file

#Read the forecast into a weather variable
read weather < $file

#Twit the weather variable away
curl --basic --user $twitbot --data status="$weather" http://twitter.com/statuses/update.xml
```

### Batch example

Here's the [marvin\_bot](http://twitter.com/marvin_bot) agent script in batch. However, previously, you are required to have [cygwin](http://cygwin.com/) installed and the `%CYGWIN%\bin` path added to your `%PATH%` environment variable.

```
@ECHO off
REM MARVIN BOT script

REM Goto here
set here="D:\Users\Guillaume\guillaume\biscates\bots"
cd %here%

REM Marvin feed address
set feed="marvinquotes_srtd.txt"

REM temporary file
set file="marvin_bot"

REM Marvin twitter bot
set user=marvin_bot
set pass=****

REM Read the feed and randomize it
cat %feed% | shuf | head -n 1 > %file%.txt

REM Read the message into a weather variable
set /P mess=<%file%.txt
echo %mess% >> %file%.log

REM Timestamp the log file
echo -------------------------------------------->> %file%.log
echo %DATE% %TIME% >> %file%.log

REM Twit the weather variable away
wget --http-user=%user% --http-pass=%pass% --post-data="status=%mess%" -O %file%.log http:// twitter.com/statuses/update.xml
@ECHO on
```

### Check if message is new

This code in bash or batch checks for new messages to be twittered.

#### bash

```
touch "mess.log"
if [ `diff mess.tmp mess.log | grep -e "^<" | wc -l` -eq 1 ]
then

 #Uncomment one the following couple of lines
 #cat mess.tmp >> mess.log #keeps historical record of all messages
 #cat mess.tmp > mess.log #keeps only last message

 #...

fi
```

#### batch

NOTE: you need to find windows Xp [binaries](http://unxutils.sourceforge.net/) of grep, diff, sed, (g)awk, tail, head, cat, tac, sort and install them in your PATH.

```
touch "mess.old"
diff mess.tmp mess.old | grep -e "^<" | wc -l > diff.tmp
set /P diff=<diff.tmp
if %diff% neq 1 goto end
goto mess

:mess
REM >> keeps historical record
REM > keeps last mess
REM Uncomment one the following couple of lines
REM type mess.tmp>>mess.old
REM type mess.tmp>mess.old

:end
```

### Sample food

```
> cat Novembro_SAASUTL.htm | \
perl -ne 'm/>([^<>]*?)<\// && print$1."\n"' | \
sed -f sasutl.sed -e '/^$/d' | \
awk '/29-11-2007/,/30-11-2007/'
```

1. Get [http://www.sas.utl.pt/alim/Ementas/View/View.aspx](http://www.sas.utl.pt/alim/Ementas/View/View.aspx) to \\\\guillaume\\temp
2. Transfer via winscp to point.pt -p 2996

Status update

## External References

- [Twitter home](http://www.twitter.com)
- [Twitter api](http://groups.google.com/group/twitter-development-talk/web/api-documentation)
- [GNU for windows](http://unxutils.sourceforge.net/).
