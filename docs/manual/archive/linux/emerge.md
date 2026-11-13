# Emerge

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Emerge is a package manager typically found on Gentoo. Use it to install missing programs and required libraries.

## Searching for packages

```
>emerge --search string
```

## Installing a package

```
>emerge package
```

## Example

```
>emerge --search rpm
Searching...
[ Results for search key : rpm ]
[ Applications found : 4 ]

- app-arch/deltarpm [ Masked ]
      Latest version available: 3.3
      Latest version installed: [ Not Installed ]
      Size of downloaded files: 460 kB
      Homepage:    ftp://ftp.suse.com/pub/projects/deltarpm/
      Description: tools to create and apply deltarpms
      License:     BSD

- app-arch/rpm
      Latest version available: 4.2-r1
      Latest version installed: 4.2-r1
      Size of downloaded files: 28,541 kB
      Homepage:    http://www.rpm.org/
      Description: Red Hat Package Management Utils
      License:     GPL-2 LGPL-2

- app-arch/rpm2targz
      Latest version available: 9.0-r3
      Latest version installed: 9.0-r3
      Size of downloaded files: 7 kB
      Homepage:    http://www.slackware.com/config/packages.php
      Description: Convert a .rpm file to a .tar.gz archive
      License:     as-is

- dev-perl/RPM
      Latest version available: 0.40-r1
      Latest version installed: [ Not Installed ]
      Size of downloaded files: 53 kB
      Homepage:    http://search.cpan.org/author/RJRAY/Perl-RPM-0.40/
      Description: RPM:: module for perl
      License:     || ( Artistic GPL-2 )
```

```
>emerge app-arch/rpm
```
