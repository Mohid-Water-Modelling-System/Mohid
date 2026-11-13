# CopyMohidProject

Eventually it is useful to copy a full aplication of Mohid but without the output (the hdf5 files).

## Copy (win)

To do so is fairly easy from the command line and with the addition of an exclude file types file. Here's what you should type:

```
xcopy /C /E .\CurrAplication\*.* .\NewAplication /EXCLUDE:exclude.txt
```

That command will create a duplicata of folder *CurrAplication* but without the files whose extensions are cited in *exclude.txt*. The exclude.txt file should contain the following text:

```
.hdf5
.srh
.bmp
.fin
.srw
.ppt
.jpg
.gif
.pdf
.ps
.eps
.doc
.log
.zip
```

## Configure (win)

Now you need to configure your project so it can run from Mohid GUI. To do so you need to edit the *.moh* file and replace all the references to *CurrAplication* with *NewAplication.*
