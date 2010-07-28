using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;

//using Common.Logging;


namespace MOHID.OpenMI.MohidLand.Wrapper.UnitTest
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {

            System.Collections.Hashtable ht = new System.Collections.Hashtable();
            ht.Add("FilePath", @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Catchment\exe\nomfich.dat");
            MohidLandEngineWrapper w = new MohidLandEngineWrapper();
            w.Initialize(ht);

            //MohidLandEngineDotNetAccessTest test = new MohidLandEngineDotNetAccessTest();
            //test.Init();
            //test.GetModelID();
            //test.AccessTimes();
            //test.RunWholeSimulation();
            //test.ClearUp();
        }
    }
}
