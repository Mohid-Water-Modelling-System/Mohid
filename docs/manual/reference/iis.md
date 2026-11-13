# IIS

**IIS** is the windows web server service that runs web-pages.

## How to execute shell scripts

Should the web-page refuse to execute some command-line script, then it is likely that the "cmd.exe" doesn't has permissions for the IIS user session to access it. Here's how to fix it:

1. Right-click "cmd.exe" from the windows folder on Properties.
2. Select *security* --> *add* --> *advanced* --> *locations* --> choose the computer --> type 'IUS' --> *check names* --> select "IUSR\_MACHINE" --> Enable Read/Execute permissions --> Apply/Ok --> close.
3. Restart the IIS service.
