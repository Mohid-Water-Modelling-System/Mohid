# Scp

Secure Copy for linux systems,

```
>scp mpich.tar.gz guillaume@nfist.ist.utl.pt:mpich.tar.gz
```

This above line copies a local file to a target host.

```
>scp guillaume@nfist.ist.utl.pt:mpich.tar.gz fedora@localhost
```

The above line does the converse from the prior.

## How to use scp without prompting for passwords

NOTE: CITING FROM [THIS](http://www.linuxjournal.com/article/8600) REFERENCE

Here's an example demonstrating the most basic syntax for scp. To copy a file named abc.tgz from your local PC to the /tmp dir of a remote PC called bozo, use:

```
scp abc.tgz root@bozo:/tmp
```

You now are asked for bozo's root password, so we're not quite there yet. The system still is asking for a password, so it's not easily scriptable. To fix that, follow this one-time procedure, after which you can make endless password-less scp copies:

Decide which user on the local machine will be using scp later on. Of course, root gives you the most power, and that's how I personally have done it. I'm not going to give you a lecture here on the dangers of root, so if you don't understand them, choose a different user. Whatever you choose, log in as that user now and stay there for the rest of the procedure. Log in as this same user when you use scp later on.

Generate a public/private key pair on the local machine. Say what? If you're not familiar with public key cryptography, here's the 15-second explanation. In public key cryptography, you generate a pair of mathematically related keys, one public and one private. You then give your public key to anyone and everyone in the world, but you never ever give out your private key. The magic is in the mathematical makeup of the keys; anyone with your public key can use it to encrypt a message, but only you can decrypt it with your private key. Anyway, the syntax to create the key pair is:

```
ssh-keygen -t rsa
```

In response, you should see:

```
Generating public/private rsa key pair
Enter file in which to save the key ...
```

```
Press Enter to accept this.
```

In response, you should see:

```
Enter passphrase (empty for no passphrase):
```

You don't need a passphrase, so press Enter twice.

In response, you should see:

```
Your identification has been saved in ...
Your public key has been saved in ...
```

Note the name and location of the public key just generated. It always ends in .pub.

Copy the public key just generated to all of your remote Linux boxes. You can use scp or FTP or whatever to make the copy. Assuming you're using root--again, see my warning in step 1--the key must be contained in the file /root/.ssh/authorized\_keys. Or, if you are logging in as a user, for example, clyde, it would be in /home/clyde/authorized\_keys. Notice that the authorized\_keys file can contain keys from other PCs. So, if the file already exists and contains text, you need to append the contents of your public key file to what already is there.

Now, with a little luck, you should be able to scp a file to the remote box without needing to use a password. So let's test it by trying our first example again. Copy a file named xyz.tgz from your local PC to the /tmp dir of a remote PC called bozo:

```
scp xyz.tgz root@bozo:/tmp
```

Wow--it copied with no password!

A word about security before we go on. This local PC just became pretty powerful, as it now has access to all the remote PCs with only the one local password. So that one password better be strong and well guarded.
