# Alias

**Alias** command in Linux allows to define customized commands or shortcuts in linux shell sessions.

## How to make alias commands permanent

**Alias** commands can be made permanent for a user's session. Simply edit your *.bashrc* (for bash sessions) file in your home directory and add all your alias commands. Save the file and restart a new terminal session. That's it.

Here's a snippet from a *[.bashrc](bash.md)* file

```
# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.

# some more ls aliases
alias ll='ls -l -a'
alias la='ls -A'
alias l='ls -CF'
alias maretec='cd /home/user/certifs; sudo openvpn --config maretec_ist.ovpn;'
```

## Example

```
> alias ll='ls -l -a'
```

Will create the `ll` shortcut for the rest of the session.

```
> alias
```

will list the user's available aliases.

## External references

- [Permanent aliases](http://www.bellevuelinux.org/make_alias_permanent.html)
