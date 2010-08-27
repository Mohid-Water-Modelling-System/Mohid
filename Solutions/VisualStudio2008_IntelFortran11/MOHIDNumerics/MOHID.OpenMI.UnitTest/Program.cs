using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using OpenMI.Standard;
using Oatc.OpenMI.Sdk.Backbone;
using MOHID.OpenMI.MohidWater.Wrapper;
using MOHID.OpenMI.MohidLand.Wrapper;
using Oatc.OpenMI.Sdk.DevelopmentSupport;

//using Common.Logging;


namespace MOHID.OpenMI.UnitTest
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            //Runs a test of Mohid Land (usefull for debug)
            runMohidLand();

            //Runs a test of Mohid Water(usefull for debug)
            runMohidWater();

            //System.Collections.Hashtable ht = new System.Collections.Hashtable();
            //ht.Add("FilePath", @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Catchment\exe\nomfich.dat");
            //MohidLandEngineWrapper w = new MohidLandEngineWrapper();
            //w.Initialize(ht);
            //w.Finish();
            //w.Dispose();
            //w.Initialize(ht);
            //w.Finish();
            //w.Dispose();


            //MohidLandEngineDotNetAccessTest mohidLandTest = new MohidLandEngineDotNetAccessTest();

            //mohidLandTest.Init();
            ////mohidLandTest.RunWholeSimulation();
            //mohidLandTest.ClearUp();
            //mohidLandTest.Init();
            //mohidLandTest.ClearUp();


            //MohidWaterEngineDotNetAccessTest mohidWaterTest = new MohidWaterEngineDotNetAccessTest();

            //mohidWaterTest.Init();
            ////mohidWaterTest.GetModelID();
            //mohidWaterTest.ClearUp();

            //mohidWaterTest.Init();
            //mohidWaterTest.ClearUp();


        }

        private static void runMohidLand()
        {
            System.Collections.Hashtable ht = new System.Collections.Hashtable();
            ht.Add("FilePath", @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Catchment\exe\nomfich.dat");
            MohidLandEngineWrapper w = new MohidLandEngineWrapper();
            w.Initialize(ht);

            ITimeSpan modelSpan = w.GetTimeHorizon();
            double now = modelSpan.Start.ModifiedJulianDay;

            double downStreamWaterLevel = 0.0;

            
            while (now < modelSpan.End.ModifiedJulianDay)
            {

                DateTime currentTime = CalendarConverter.ModifiedJulian2Gregorian(now);
                DateTime intitalTime = CalendarConverter.ModifiedJulian2Gregorian(w.GetTimeHorizon().Start.ModifiedJulianDay);

                downStreamWaterLevel = (currentTime - intitalTime).TotalHours/10.0;

                Console.WriteLine(currentTime.ToString());
                w.PerformTimeStep();

                //Gets outputs Items
                for (int i = 0; i < w.GetOutputExchangeItemCount(); i++)
                {
                    OutputExchangeItem ouputItem = w.GetOutputExchangeItem(i);

                    IValueSet values = w.GetValues(ouputItem.Quantity.ID, ouputItem.ElementSet.ID);

                }

                //Sets Input Items
                for (int i = 0; i < w.GetInputExchangeItemCount(); i++)
                {
                    InputExchangeItem inputItem = w.GetInputExchangeItem(i);

                    double[] aux = new double[1];
                    aux[0] = downStreamWaterLevel;
                    IValueSet values = new ScalarSet(aux);

                    w.SetValues(inputItem.Quantity.ID, inputItem.ElementSet.ID, values);

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

            double flow = 0.0;

            while (now < modelSpan.End.ModifiedJulianDay)
            {

                flow = flow + 0.1;

                //Gets Output exchange items
                for (int i = 0; i < w.GetOutputExchangeItemCount(); i++)
                {
                    OutputExchangeItem outputItem = w.GetOutputExchangeItem(i);

                    IValueSet values = w.GetValues(outputItem.Quantity.ID, outputItem.ElementSet.ID);
                }

                //Sets Input Items
                for (int i = 0; i < w.GetInputExchangeItemCount(); i++)
                {
                    InputExchangeItem inputItem = w.GetInputExchangeItem(i);

                    double[] aux = new double[1];
                    aux[0] = flow;
                    IValueSet values = new ScalarSet(aux);

                    w.SetValues(inputItem.Quantity.ID, inputItem.ElementSet.ID, values);
                   
                }


                w.PerformTimeStep();


                now = w.GetEarliestNeededTime().ModifiedJulianDay;

            }


            w.Finish();
        }
    }
}
