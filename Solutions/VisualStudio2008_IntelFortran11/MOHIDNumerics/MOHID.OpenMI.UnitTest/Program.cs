using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using OpenMI.Standard;
using Oatc.OpenMI.Sdk.Backbone;
using MOHID.OpenMI.MohidWater.Wrapper;

//using Common.Logging;


namespace MOHID.OpenMI.MohidLand.Wrapper.UnitTest
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            //Runs a test of Mohid Land (usefull for debug)
            //runMohidLand();

            //Runs a test of Mohid Water(usefull for debug)
            runMohidWater();

        }

        private static void runMohidLand()
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

        private static void runMohidWater()
        {
            System.Collections.Hashtable ht = new System.Collections.Hashtable();
            ht.Add("FilePath", @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Estuary\exe\nomfich.dat");
            MohidWaterEngineWrapper w = new MohidWaterEngineWrapper();
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

                    //IValueSet values = w.GetValues(ouputItem.Quantity.ID, ouputItem.ElementSet.ID);

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
