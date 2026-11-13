# TortoiseSvn

**TortoiseSvn** is a very cool [svn](../archive/linux/svn.md) client for windows with a cool interface that blends seamlessly with the windows explorer.

## How to use the svn client with Mohid

Simply install presto the tortoise svn [Tortoise svn](http://tortoisesvn.tigris.org/). You can also see [video screencasts](../archive/web/codeplex.md) instead of the picture-only tutorial below...

### Creating a new local repo

Here's the quick how-to add a new **mohid\_net** project from our local intranet repository:

- Create a new folder:

- Right click the new folder and select **svn checkout...**:

- Set the central server repo address as **[https://mohid.svn.codeplex.com/svn](https://mohid.svn.codeplex.com/svn)**, then click **ok**:

- Insert your subversion user name and password. Then click **ok**.

- Get the files, and when it's done click **ok**:

- That's it! You've imported the whole project!

### Getting the latest version of the code

- To refresh and get the latest version of the code simply right-click on your local-repository, and then select the **SVN** update option:

- Once the new files are downloaded, simply hit the **ok** button.

- That's it!

### Commiting changes to the code

- Now whenever you feel like to commit your changes, if any (provided that you are granted the rights to do so - ask the mohid code administrators for a user name and the write permission), first you need to update and merge to the last version of the code. Simply right click the subversioned root folder and click **Svn update...**. Get the latest version and then beware if any merging occurs. Usually, 90% of the time, the merges made are automatic. Sometimes however you may be prompted by the program to merge manually the changed line of code.

- After all the merges are made, please make sure that your modifications have prevailed in your local copy of the code and that your project still compiles successfully.

- After the merge is made and the compiling works successfully, you can finally proceed to commit the code: simply right click the subversioned folder, and click **Svn commit...**:

- Make sure you log all your changes, then click **ok**:

- And... it's done!

## External References

- <goto>TortoiseSvn</goto>
- <goto>RapidSvn</goto>
