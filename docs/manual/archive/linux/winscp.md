# Winscp

**Winscp** is a little yet powerful secure copy ([scp](scp.md)) or secure file protocol transfer (sftp) client. In command line mode it allows to easily automate some tasks such as access and transfer from remote linux-based machines.

## Examples

The following line executes a text-script that transfers to the machine *192.168.20.160* the folder *MyFolder*:

```
>winscp3 /console /script=winscpscript.txt
```

```
#winscript.txt content
option batch on
option confirm off
open fedora:fedora@192.168.20.160
cd /home/fedora/Aplica
option transfer binary
put MyFolder
close
exit
```

## Applications

Here's a combination with a batch file, that creates a clean copy of a MOHID simulation and then transfers it to a linux machine with a winscp script:

```
> copysimulation MySimulation
```

```
#copysimulation.bat content
#Input argument %1 stands for simulation name
@echo off
set COPYHOME=C:\Program Files\2gif
xcopy %1\*.* %1_copy /E /EXCLUDE:%COPYHOME%\exclui.txt
perl %COPYHOME%/replace.pl 's/%2/\.\.\/\.\./g' *
perl %COPYHOME%/replace.pl 's/\\/\//g' *
winscp3 /console /script=%COPYHOME%\winscpfedora.txt
@echo on
```

```
#winscpfedora.txt content
#don't stall on errors
option batch on
#don't ask for confirmation
option confirm off

#connect
open fedora:fedora@192.168.23.151
cd /homes/Aplica
option transfer binary
put

#disconnect
close
exit
```

## External references

- The [download page](http://winscp.net/eng/download.php#download2).
- The documented [winscp site](http://winscp.net/eng/docs/scripting).
