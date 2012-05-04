using System;


namespace MOHID.OpenMI.UnitTest
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            //Runs a test of Mohid Land (usefull for debug)
            //runMohidLand();

            //Runs a test of Mohid Water(usefull for debug)
            //runMohidWater();

            ////Runs Integrates Tests
            //IntegratedTests integratedTests = new IntegratedTests();
            //integratedTests.Test_T2();

        }

        private static void runMohidLand()
        {

            MohidLandEngineTests mohidLandEngineTests = new MohidLandEngineTests();
            mohidLandEngineTests.Init();
            mohidLandEngineTests.RunSimulationWithInputAndOutput();
            mohidLandEngineTests.ClearUp();

            //mohidLandEngineTests.Init();
            //mohidLandEngineTests.GetModelID();
            //mohidLandEngineTests.ClearUp();

            //mohidLandEngineTests.Init();
            //mohidLandEngineTests.AccessTimes();
            //mohidLandEngineTests.ClearUp();

  
        }

        private static void runMohidWater()
        {
            MohidWaterEngineTests mohidWaterEngineTests = new MohidWaterEngineTests();
            mohidWaterEngineTests.Init();
            mohidWaterEngineTests.RunSimulationWithInputAndOutput();
            mohidWaterEngineTests.ClearUp();

            //mohidWaterEngineTests.Init();
            //mohidWaterEngineTests.GetModelID();
            //mohidWaterEngineTests.ClearUp();

            mohidWaterEngineTests.Init();
            mohidWaterEngineTests.AccessTimes();
            mohidWaterEngineTests.ClearUp();

        }



    }
}
