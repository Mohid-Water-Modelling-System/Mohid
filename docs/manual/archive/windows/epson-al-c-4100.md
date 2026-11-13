# Epson AL C-4100

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

The EPSON is a fine printer, but installing the drivers in Win XP 64 can be troublesome.

## Installing drivers

1. Get the drivers *\\\\Davinci\\Software\\ImpressoraEpson\_Drivers* and copy the folder to somewhere in your computer (e.g. Software\\ImpressoraEPSON\_Drivers),
2. From the control panel do an Add printer,
3. Select *Local Printer Attached to This Computer*. Uncheck *automatically detect my printer's drivers*,
4. *Create a New Port*\-->*Standard TCP/IP*,
5. Printer IP --> 192.168.20.128,
6. *Have disk ...* and select the .INF file in the folder in your computer where you saved the drivers (e.g. Software\\ImpressoraEPSON\_Drivers),
7. Say *ok* to everything,
8. That's it!

## Troubleshooting

### No connection

The other day, we lost connection with the printer although she seemed to be working fine. The problem was its ip, which usually was on 192.168.20.127, somehow, changed to 192.168.20.169. And its network mode was set in AUTO. Filipa, brought it back to work by printing *Network Status Sheets* from the printer's menu and by changing the network mode to PANEL and then setting again its old ip, 192.168.20.128, again, from the printer's menu.

More Details:

1. Go to Network menu->Network Config
2. Change Network Config from "No" to "Yes".
3. Change Get IP Address from "Auto" to "Panel".
4. Change the IP with the enter (to travel between numbers) and with up and down (to change values). Put the default 192.168.20.127
5. After this go to Reset Menu and do a "Reset all"

## External References

[EPSON support](http://esupport.epson-europe.com/ProductHome.aspx?lng=pt-PT)
