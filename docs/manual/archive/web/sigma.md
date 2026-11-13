# Sigma

**Sigma** cluster is where personal web-pages and group web-pages of ist.utl.pt are stored.

## Quick login

### Windows

1. Download [putty](http://the.earth.li/~sgtatham/putty/latest/x86/putty.exe),
2. Fill in the **Host** field with **sigma.ist.utl.pt**, then click on open,
3. Insert your fenix login without the *@ist.utl.pt* part, then your password,

You're in!

### Linux

```
> ssh sigma.ist.utl.pt -l istxxxxxx
```

## Personal web-page

```
http://web.ist.utl.pt/istxxxxxx
> cd /afs/.ist.utl.pt/users/istxxxxxx
```

## Groups web-page

```
http://(group).ist.utl.pt
> cd /afs/.ist.utl.pt/groups/(group)/web
```

### Example

```
http://groups.ist.utl.pt/maretec
> cd /afs/.ist.utl.pt/groups/maretec/web
```

### Administrators

Administrators may add new group administrators by using the [ciist self-service for groups management web-page](https://ciist.ist.utl.pt/servicos/self_service/show_groups.php). The personal IST ID login is required.

## Mounting AFS network drives

You can access to your sigma AFS files through a network drive mounted on your computer. This is a more convenient way to copy/transfer files from your local machine to the AFS server.

### Windows

- [Ciist wiki for windows AFS clients](https://ciist.ist.utl.pt/ciwiki/AfsClienteWindows)

### Linux

- [Ciist wiki for linux AFS clients](https://ciist.ist.utl.pt/ciwiki/AfsClienteLinux)

## External links

- [ciist serviços](https://ciist.ist.utl.pt/servicos/web.php)
- [Maretec](http://maretec.ist.utl.pt/)
- [Maretec](http://groups.ist.utl.pt/maretec)
