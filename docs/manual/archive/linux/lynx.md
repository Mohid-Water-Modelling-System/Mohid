# Lynx

[Lynx](http://en.wikipedia.org/wiki/Lynx_%28web_browser%29) is a text-based web browser for [linux](nix-platforms.md). Useful for command line web-browsing.

This command records the keystrokes used to download a file to the current directory.

```
>lynx -cmd_log:/etc/lynx/lynx_getfile.log http://www.mohid.com/hydrogroup/presentations/file.ppt
```

This command exectuts the scripted keystrokes and downloads a file to the current directory.

```
>lynx -cmd_script:/etc/lynx/lynx_getfile.log http://www.mohid.com/hydrogroup/presentations/file.ppt
```
