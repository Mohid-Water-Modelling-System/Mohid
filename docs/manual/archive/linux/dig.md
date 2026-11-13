# Dig

**dig** requests to a dns server it's list of ips mapped to urls. The dns server asks the [bind](bind.md) service for it to indicate the correct file zone where the requested urls are located. Here's how to [configure bind](bind.md) and here's how to configure a dns zone.

## Examples

```
> dig @72.14.235.9 any
```

```
> dig www.google.com any
```

In this case, dig will first query my ip provider's dns server for the ip associated to the url *"www.google.com"* and only after it will dig for the ip's mapped by the *"www.google.com"'*s dns server.

```
> dig -x "ip"
```
