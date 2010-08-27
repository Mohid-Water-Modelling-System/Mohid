using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using MOHID.OpenMI.MohidWater.Wrapper;

namespace MOHID.OpenMI.UnitTest
{
    [TestFixture]
    public class MohidWaterEngineDotNetAccessTest
    {

        MohidWaterEngineDotNetAccess mohidWaterEngineDotNetAccess;
        String _filePath;
 
        [SetUp]
        public void Init()
        {
            mohidWaterEngineDotNetAccess = new MohidWaterEngineDotNetAccess();
            _filePath = @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Estuary\exe\nomfich.dat";
            mohidWaterEngineDotNetAccess.Initialize(_filePath);
        }
 
        [TearDown]
        public void ClearUp()
        {
            mohidWaterEngineDotNetAccess.Finish();
        }
 
        [Test]
        public void GetModelID()
        {
            try
            {
                String modelID = mohidWaterEngineDotNetAccess.GetModelID();
                Assert.AreEqual("Sample Estuary", modelID);
            }
            catch(System.Exception e)
            {
                Console.WriteLine(e.Message);
                //this.WriteException(e.Message);
                throw(e);
            }
        }


        [Test]
        public void RunWholeSimulation()
        {
            mohidWaterEngineDotNetAccess.RunSimulation();
        }

        [Test]
        public void AccessTimes()
        {
            DateTime start = new DateTime(2002, 1, 1, 0, 0, 0);
            DateTime end = new DateTime(2002, 1, 1, 12, 0, 0);

            DateTime startInstant = mohidWaterEngineDotNetAccess.GetStartInstant();
            DateTime endInstant = mohidWaterEngineDotNetAccess.GetStopInstant();
            Double timeStep = mohidWaterEngineDotNetAccess.GetCurrentTimeStep();

            Assert.AreEqual(startInstant.Ticks, start.Ticks);
            Assert.AreEqual(endInstant.Ticks, end.Ticks);
            Assert.AreEqual(15.0, timeStep);
        }
    
    }
}
