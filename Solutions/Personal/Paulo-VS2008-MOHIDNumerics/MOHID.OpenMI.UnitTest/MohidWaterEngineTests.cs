using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using NUnit.Framework;
using MOHID.OpenMI.MohidWater.Wrapper;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.DevelopmentSupport;
using OpenMI.Standard;

namespace MOHID.OpenMI.UnitTest
{
    [TestFixture]
    public class MohidWaterEngineTests
    {

        private MohidWaterEngineWrapper mohidWaterEngineWrapper;
 
        [SetUp]
        public void Init()
        {
            mohidWaterEngineWrapper = new MohidWaterEngineWrapper();
            System.Collections.Hashtable ht = new System.Collections.Hashtable();
            ht.Add("FilePath", @"D:\MohidProjects\Studio\20_OpenMI\Sample Estuary\exe\nomfich.dat");
            mohidWaterEngineWrapper.Initialize(ht);
        }
 
        [TearDown]
        public void ClearUp()
        {
            mohidWaterEngineWrapper.Finish();
            mohidWaterEngineWrapper.Dispose();
        }
 
        [Test]
        public void GetModelID()
        {
            String modelID = mohidWaterEngineWrapper.GetModelID();
            Assert.AreEqual("Sample Estuary", modelID);
        }

        [Test]
        public void AccessTimes()
        {
            DateTime start = new DateTime(2002, 1, 1, 0, 0, 0);
            DateTime end = new DateTime(2002, 1, 1, 12, 0, 0);

            ITimeSpan timeHorizon = mohidWaterEngineWrapper.GetTimeHorizon();

            Assert.AreEqual(CalendarConverter.ModifiedJulian2Gregorian(timeHorizon.Start.ModifiedJulianDay), start);
            Assert.AreEqual(CalendarConverter.ModifiedJulian2Gregorian(timeHorizon.End.ModifiedJulianDay), end);
        }

        [Test]
        public void RunSimulationWithInputAndOutput()
        {
            ITimeSpan modelSpan = mohidWaterEngineWrapper.GetTimeHorizon();
            double now = modelSpan.Start.ModifiedJulianDay;

            Stopwatch win = new Stopwatch();
            Stopwatch wout = new Stopwatch();
            Stopwatch wengine = new Stopwatch();

            while (now < modelSpan.End.ModifiedJulianDay)
            {

                DateTime currentTime = CalendarConverter.ModifiedJulian2Gregorian(now);
                DateTime intitalTime =
                    CalendarConverter.ModifiedJulian2Gregorian(
                        mohidWaterEngineWrapper.GetTimeHorizon().Start.ModifiedJulianDay);

                Console.WriteLine(currentTime.ToString());

                wengine.Start();
                mohidWaterEngineWrapper.PerformTimeStep();
                wengine.Stop();

                wout.Start();
                //Gets outputs Items
                for (int i = 0; i < mohidWaterEngineWrapper.GetOutputExchangeItemCount(); i++)
                {
                    OutputExchangeItem ouputItem = mohidWaterEngineWrapper.GetOutputExchangeItem(i);

                    IValueSet values = mohidWaterEngineWrapper.GetValues(ouputItem.Quantity.ID, ouputItem.ElementSet.ID);

                }
                wout.Stop();

                //Sets Input Items
                win.Start();
                for (int i = 0; i < mohidWaterEngineWrapper.GetInputExchangeItemCount(); i++)
                {
                    InputExchangeItem inputItem = mohidWaterEngineWrapper.GetInputExchangeItem(i);

                    double[] aux = new double[inputItem.ElementSet.ElementCount];
                    for (int j = 0; j < inputItem.ElementSet.ElementCount; j++)
                    {
                        aux[j] = 0;
                    }
                    IValueSet values = new ScalarSet(aux);

                    //mohidLandEngineWrapper.SetValues(inputItem.Quantity.ID, inputItem.ElementSet.ID, values);

                }
                win.Stop();

                now = mohidWaterEngineWrapper.GetEarliestNeededTime().ModifiedJulianDay;
            }

            Console.WriteLine("Input Exchange:  " + win.ElapsedMilliseconds.ToString());
            Console.WriteLine("Output Exchange: " + wout.ElapsedMilliseconds.ToString());
            Console.WriteLine("Engine:          " + wengine.ElapsedMilliseconds.ToString());
        }

    
    }
}
