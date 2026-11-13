# Adduser

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Adduser allows to add new users in linux. Here's the help page:

```
[root@dual_64 home]# adduser --hel
Usage: useradd [options] LOGIN

Options:
  -b, --base-dir BASE_DIR       base directory for the new user account
                                home directory
  -c, --comment COMMENT         set the GECOS field for the new user account
  -d, --home-dir HOME_DIR       home directory for the new user account
  -D, --defaults                print or save modified default useradd
                                configuration
  -e, --expiredate EXPIRE_DATE  set account expiration date to EXPIRE_DATE
  -f, --inactive INACTIVE       set password inactive after expiration
                                to INACTIVE
  -g, --gid GROUP               force use GROUP for the new user account
  -G, --groups GROUPS           list of supplementary groups for the new
                                user account
  -h, --help                    display this help message and exit
  -k, --skel SKEL_DIR           specify an alternative skel directory
  -K, --key KEY=VALUE           overrides /etc/login.defs defaults
  -m, --create-home             create home directory for the new user
                                account
  -l,                       do not add user to lastlog database file
  -M,                       do not create user's home directory(overrides /etc/login.defs)
  -r,                       create system account
   -o, --non-unique              allow create user with duplicate
                                (non-unique) UID
  -p, --password PASSWORD       use encrypted password for the new user
                                account
  -s, --shell SHELL             the login shell for the new user account
  -u, --uid UID                 force use the UID for the new user account
```

## Adding users in Dual\_64

```
> adduser -m USERNAME
> passwd USERNAME
passwd> PASSWORD
> chmod -R --reference=/home/fedora /home/USERNAME
> ln -s /home/Aplica /home/USERNAME
> ln -s /home/Projects /home/USERNAME
> ln -s /home/Software /home/USERNAME
```

## Removing users

```
> userdel -r USERNAME
```

## Adding sudoers

Edit the */etc/sudoers*

```
root> visudo
visudo> USERNAME    ALL=(ALL)    ALL
```
