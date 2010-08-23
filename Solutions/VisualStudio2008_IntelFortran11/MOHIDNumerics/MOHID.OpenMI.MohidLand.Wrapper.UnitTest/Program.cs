using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using OpenMI.Standard;
using Oatc.OpenMI.Sdk.Backbone;

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

            ITimeSpan modelSpan = w.GetTimeHorizon();
            double now = modelSpan.Start.ModifiedJulianDay;
            while (now < modelSpan.End.ModifiedJulianDay)
            {
                w.PerformTimeStep();

                //Gets outputs Items
                for (int i = 0; i < w.GetOutputExchangeItemCount(); i++)
                {
                    OutputExchangeItem ouputItem = w.GetOutputExchangeItem(i);

                    IValueSet values = w.GetValues(ouputItem.Quantity.ID, ouputItem.ElementSet.ID);

                    //if (values is ScalarSet)
                    //{
                    //    Console.WriteLine(((ScalarSet)values).data[0].ToString());
                    //}

                }
                
                now = w.GetEarliestNeededTime().ModifiedJulianDay;
            }


            w.Finish();


        }
    }
}
