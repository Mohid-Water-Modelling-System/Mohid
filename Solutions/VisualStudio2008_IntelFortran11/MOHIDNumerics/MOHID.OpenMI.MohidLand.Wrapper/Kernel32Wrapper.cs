using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace MOHID.OpenMI.MohidLand.Wrapper
{
    /// <summary>
    /// Summary description for Kernel32Wrapper.
    /// </summary>
    public class Kernel32Wrapper
    {

        [DllImport("kernel32.dll")]
        public static extern bool FreeLibrary(IntPtr hModule);

        [DllImport("kernel32.dll")]
        public static extern IntPtr LoadLibrary(string lpFileName);

    }
}
