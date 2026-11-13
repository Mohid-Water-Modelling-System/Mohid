# Timestamp ISO conventions

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

A Timestamp convention could be *13th August 2006* or *13/08/2006*, or yet *08/13/2006* followed by the hours: *13:30.5* or *13:30:30* or *13h03m30s* etc...

## Recommended Timestamp Convention

For example if one would choose the *13-08-2006* convention then here's a date set's chronological sorting:

- 13-08-2006
- 11-09-2006
- 12-10-2006

And here's the same set's alphabetical sorting (this is the default sorting when prompting for a directory listing):

- 11-09-2006
- 12-10-2006
- 13-08-2006

The alphabetical and chronological sorting don't match! This is a nuisance, especially when dealing with timestamped files (ex: *11-09-2006\_myfile.dat*). Thus, the best suited Timestamp Convention in computerized systems is one that matches alphabetical sorting with chronological sorting. The simplest one is the *20060813* convention. Here's the same set's alphabetical sorting:

- 20060811
- 20060812
- 20060813

and here's the chronological sorting:

- 20060811
- 20060812
- 20060813

They're the same! Thus, throughout this wiki, editors are encouraged to use the *20060813* timestamp convention.
