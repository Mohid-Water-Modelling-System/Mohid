# Switching from SourceSafe to Subversion

## FAQ

### *Are we really going to change from source-safe to subversion?*

Yes, after a lot of propositions and debating with the main developers at Maretec and Hidromod, we opted to make the switch from source-safe to subversion.

### *Do I really have to change?*

Yes.

### *Why?*

Because source-safe has some limitations that prevent a better development and maintenance of the MOHID code:

- Source-safe is proprietary code (you need to pay for a licence to use it),
- Source-safe is no longer supported and has been discontinued by the maker,
- Source-safe is hardly portable to other platforms such as linux or mac.

On the other hand:

- Subversion is free-software,
- Subversion has a very active community of users and developers,
- Subversion is cross-platform (windows, linux, mac, unix, solaris...),

### *When?*

In the next coming months. June 2009 is a good deadline. We'll let you know in due time. Documentation will be made and workshops will be held.

### *What are the implications?*

The users of source-safe will have to install a subversion client in their computers and configure it in order to access the new subversion repository of the MOHID code. So, basically, it's a new interface to get the latest version of the code. And that's all.

### *What are the other implications?*

Developers will have to use the subversion client to perform all the standard source-safe actions like *check-out*, *check-in*, *diff*, *merge*, *roll-back* and *branching*. Also, the same code can be edited at the same time. The latest to check-in will have to undergo the *merge* process carefully. The experience of the subversion user's database tells that it's better that way.

### '*Do all the versions of the code will be stored in Subversion?*

No. We'll start from a clean recent version. Only the subsequent versions will stored in Subversion. However, we will keep the the source-safe database unchanged and locked in "read-only".

### *Will it work the same way as with source-safe?*

Yes. *Get latest version*, *check-in*, *check-out*. Same thing, different interface.

### **Is this the right time change?**

As good as any other.
