# Crontab

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Crontab is the tool to schedule tasks in linux

```
>crontab -e
```

Example:

```
0	4	12	3	*	aplication &
```

The above example tells to execute the application "aplication" the 12/3 at 04:00. The ending ampersand (&) is to [wikipedia:fork](http://en.wikipedia.org/wiki/fork) [fork](http://en.wikipedia.org/wiki/Fork_%28operating_system%29) the application in an external shell.

## Disable email

Add this at the end of the cronjob line:

```
>/dev/null 2>&1
```

## First friday of every month

This one's a bit more complicated but it's manageable:

```
30 13 1-7 * * [ `date +\%w` -eq 5 ]&& myscript.sh >/dev/null 2>&1
```

## External References

- [http://www.devin.com.br/eitch/crontab/](http://www.devin.com.br/eitch/crontab/)
