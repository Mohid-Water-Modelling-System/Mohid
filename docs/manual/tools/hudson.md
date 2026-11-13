# Hudson

Hudson allows automatic daily builds and daily tests with a reporting system on top. Good to use with MOHID.

## Hidromod documentation

Continuous Integration Server for Mohid with Hudson by ***Eduardo Aires*** - Hidromod

In this guide I'm going to show how to set up the Mohid project on the Continuous integration server Hudson ( [http://hudson-ci.org](http://hudson-ci.org) ).

### Initial downloads

The following file are needed besides Java (at least 1.5). Get the latest version of all files and notice that the Hudson file has the extension .war and plugins .hpi.

```
Hudson server application - http://hudson-ci.org/latest/hudson.war
```

### Installation steps

I'm going to install Hudson into c:\\Hudson.

Copy the hudson.war file to c:\\Hudson

In command line go to c:\\Hudson and start Hudson through "java -DHUDSON\_HOME=data -jar hudson.war". Verify that you can access Hudson on web browser through [http://localhost:8080](http://localhost:8080)

(To stop Hudson press Ctrl+C in the command prompt where you started Hudson)

#### To start it at boot time

Create this [batch](../building/windows/batch.md) file in the *c:\\hudson* folder *StartHudson.bat*:

```
c:
cd c:\Hudson
java -DHUDSON_HOME=data -jar hudson.war
pause
```

Then add it to your windows 7 "all users" *Startup* folder in the *Start-->All Programs-->Startup* menu (right-click on the folder and select *Open all users*).

### Mohid job configuration

Click the "New job" link at the home page.

Enter the name "Mohid", check the "Build a free-style software project" and press OK.

### Source code management

Mohid uses a Subversion SCM. Hudson supports CVS and SVN out of the box, but there are many plugins for other SCMs. After checking out the files from the repository Hudson will show the new change sets since the previous build in the Build page. A detailed view of change sets can be seen in the Changes page as the name of the developer, files that were changed and the comment for the change. Each change set is linked to the Mohid subversion repository browser, so it is easy to browse the actual file that changed.

Press the Subversion radio button to configure the SCM. Repository URL= [https://mohid.svn.codeplex.com/svn](https://mohid.svn.codeplex.com/svn)

Local module directory=.

Use update=true

```
(Revert=true, for your consideration)
```

To test the configuration, press Save and then Build. The source code will be downloaded from the repository and put into the Workspace. If there were any changes in the SCM repository they can be viewed in the Changes page.

When the build is completed verify that it has checked out the code by going to the Workspace page. Using the Workspace page you can browse and view the files that has been checked out and it doesn't matter if the files are on the master or on a distributed slave!

### Build Triggers

One option is Build periodically, that execute the build processes, periodically by cron expression. Another option id Poll SCM, that ask periodically, by cron expression, if the repository has changes.

Select Poll SCM Schedule=\* \* \* \* \* It checks every minute if it contains changes, if yes do update and build processes selected after.

There is more options that can be importante: In Advanced Project Options, Quiet period – when a change have been detected, wait this amount of seconds to start update, this is for avoid multi commits by developers problems… Retry Count – if there is an error try this amount of times before send error warnings.

### Building the project

Im going to build a Debug version of MOHIDBase1, MOHIDBase2 and MOHIDWater in the MOHIDNumerics.sln solution file. Hudson also supports MSBuild, NAnt and many other build tools through plugins.

Click the "Add build step" and select "Execute Windows batch command."

Command:

```
"c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE\devenv"  "C:\Hudson\workspace\svn\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MOHIDNumerics.sln" /rebuild  "Debug|x64" /project MOHIDBase1

"c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE\devenv"  "C:\Hudson\workspace\svn\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MOHIDNumerics.sln" /rebuild  "Debug|x64" /project MOHIDBase2

"c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE\devenv"  "C:\Hudson\workspace\svn\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MOHIDNumerics.sln" /rebuild  "Debug|x64" /project MOHIDWater
```

To test the configuration, press Save and then Build. Now the source code should be updated if there any changes and then built. While the build is running you can check the Console log that is updated as the build continues.

### Integration Tests

To check if the tests don’t changed the data by mistake, we’ve created a C# application (see attached compressed file) to compare the data of two time series. The first resulted from the execution of recently build MohidWater, over a Mohid project. And the second is the previously generated for reference.

This C# application will be described latter.

To add this application or another created by you.

Click the "Add build step" and select "Execute Windows batch command."

```
c:\Hudson\Tester\mohidTester.exe c:\Hudson\Tests\Tests.xml
```

The Tests.xml is where test will be defined.

To test the configuration, do the same as in the step before. If there is an error in integration test, you can check what is different in the console log.

### Advanced Hudson

You can install hudson as windows service in [http://localhost:8080/manage](http://localhost:8080/manage) You can install a lot of plugins in http:// localhost:8080/pluginManager/ Hudson is on full speed evolution, and almost every week it have an update, try to have it updated. It can send emails with the reports, it uses rss feeds, it can use CruiseControl.Net CCTray, an system tray checker for clients [http://serverName:8080/cc.xml](http://serverName:8080/cc.xml). You can have user accounts for control Hudson, you can see statistics and much more stuff.

### mohidTester

The executable process: 1. Check if the Test.xml file passed through argument is valid. 2. For each Level in configuration from BaseDirectory will a. Copy from that level data folder Nomfich\_”RunNumber”.dat to exe folder Nomfich.dat b. Execute MohidExecutable for that level exe folder c. Validate if needed files for compare have been created d. For each ResultCompare will compare from res folder to the respective in ReferenceBaseDirectory i. Number of parameters ii. For each parameter 1. Max value 2. Min value 3. Begin instant 4. End instant 5. Model domain name 6. Name without parameter 7. Parameter 8. Parameter name 9. Number of time series values 10. For each time series values a. Compare datetime, id and value

Tests.xml

```
<?xml version="1.0" encoding="utf-8"?>
<MohidTests xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"  xmlns:xsd="http://www.w3.org/2001/XMLSchema">
<MohidExecutable>c:\Hudson\workspace\svn\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MOHI DWater\Debug\MOHIDWater.exe</MohidExecutable>
  <Projects>
    <Project>
      <Name>Sado</Name>
      <BaseDirectory>c:\Hudson\Tests\Test1\Sado</BaseDirectory>
      <ReferenceBaseDirectory>c:\Hudson\Tests\Test1\SadoReference</ReferenceBaseDirectory>
      <Runs>
        <Run>
          <Name>run1</Name>

          <Level>
            <Name>level1</Name>
            <LevelBaseDirectory>\Sado_LowResolution</LevelBaseDirectory>
            <RunNumber>1</RunNumber>
            <ResultsCompare>
              <ResultCompare>
                <ResultTest>\Station 5.srw</ResultTest>
              </ResultCompare>
            </ResultsCompare>

          </Level>

        </Run>
      </Runs>
    </Project>
   < !--
    <Project>
    ...
    </Project>
    -->
  </Projects>
</MohidTests>
```

## Configuring SMTP and Jabber server on Hudson

Go to the main Hudson settings --> Configure:

Email

```
server: smtp.gmail.com
port:  465 or 587
authentication?: yes
   user: mailing.maretec@gmail.com
   pass: the usual
   ssl?: yes
```

Jabber

```
Jabber ID: mailing.maretec@gmail.com
pass: the usual
Advanced
   server: talk.google.com
```

### Sending commands from Gtalk to Hudson

Example:

```
!build 'My project'
```

Other commands:

```
status: displays the status of the build i.e. last build result and last succesful build
health: displays the build health as computed from the last build results, test results, (compiler  warnings, findbugs warnings, ...)
build: schedules a new build for a single job
abort: abort a currently running build
help: get a list of all available commands
```

## List of available servers

Servidor na hidromod com teste de integração para mohid water [http://matrix:8080/](http://matrix:8080/)

Um system tray que verifica como está o servidor [http://nomeDoServidor:8080/cc.xml](http://nomeDoServidor:8080/cc.xml) CruiseControl.Net CCTray

## External links

- [Hudson](http://hudson-ci.org/)
- [Help](http://redsolo.blogspot.com/2008/04/guide-to-building-net-projects-using.html)
