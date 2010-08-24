using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace MOHID.OpenMI.MohidLand.Wrapper.UnitTest
{
    [TestFixture]
    public class MohidLandEngineDotNetAccessTest
    {

        MohidLandEngineDotNetAccess mohidLandEngineDotNetAccess;
        String _filePath;
 
        [SetUp]
        public void Init()
        {
            mohidLandEngineDotNetAccess = new MohidLandEngineDotNetAccess();
            _filePath = @"D:\MohidProjects\Studio\03_MOHID OpenMI\Sample Catchment\exe\nomfich.dat";
            mohidLandEngineDotNetAccess.Initialize(_filePath);
        }
 
        [TearDown]
        public void ClearUp()
        {
            mohidLandEngineDotNetAccess.Finish();
        }
 
        [Test]
        public void GetModelID()
        {
            try
            {
                String modelID = mohidLandEngineDotNetAccess.GetModelID();
                Assert.AreEqual("MOHID Land Model", modelID);
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
            mohidLandEngineDotNetAccess.RunSimulation();
        }

        [Test]
        public void AccessTimes()
        {
            DateTime startInstant = mohidLandEngineDotNetAccess.GetStartInstant();
            DateTime endInstant = mohidLandEngineDotNetAccess.GetStopInstant();
            Double timeStep = mohidLandEngineDotNetAccess.GetCurrentTimeStep();
        }
    
    }
}
