# Mathematica

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Linux

Here's how to run in an external process a mathematica batch file in \*nix systems. First, make sure you *save* the mathematica notebook file as *special-->text*.

```
> nohup nice math -noprompt -run "<<filecommand.m" > output.out &
```

You'll receive a process id. You can check the processes that you ran by using the command

```
> top
```
