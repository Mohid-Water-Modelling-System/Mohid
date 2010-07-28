using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;

namespace MOHID.OpenMI.MohidLand.Wrapper
{

    public class MohidLandEngineDotNetAccess
    {
        IntPtr _FortranDllHandle;

        /// <summary>
        /// Loads the Fortran Dll and initializes the model (calling the constructor)
        /// </summary>
        /// <param name="filePath">Path to the nomfich.dat file</param>
        public void Initialize(string filePath)
        {
            //Loads the library
            _FortranDllHandle =  Kernel32Wrapper.LoadLibrary(@"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidLandEngine\Debug\MohidLandEngine.dll");

            //Calls the constructor and reads data files, etc
            if (!(MohidLandEngineDLLAccess.Initialize(filePath, ((uint)filePath.Length))))
            {
                CreateAndThrowException();
            }
        }

        /// <summary>
        /// Performs a single time step
        /// </summary>
        public void PerformTimeStep()
        {
            if (!(MohidLandEngineDLLAccess.PerformTimeStep()))
            {
                CreateAndThrowException();
            }
        }

        /// <summary>
        /// Calls the model destructor (closes data files)
        /// </summary>
        public void Finish()
        {
            if (!(MohidLandEngineDLLAccess.Finish()))
            {
                CreateAndThrowException();
            }
            while (Kernel32Wrapper.FreeLibrary(_FortranDllHandle)) ;
        }

        public void Dispose()
        {
            if (!MohidLandEngineDLLAccess.Dispose())
            {
                CreateAndThrowException();
            }
        }

        /// <summary>
        /// Runs the whole working cycle once - Testing Only
        /// </summary>
        public void RunSimulation()
        {
            if (!(MohidLandEngineDLLAccess.RunSimulation()))
            {
                CreateAndThrowException();
            }
        }

        private void CreateAndThrowException()
        {
            int numberOfMessages = 0;
            numberOfMessages = MohidLandEngineDLLAccess.GetNumberOfMessages();
            string message = "Error Messages from MOHID Land Engine";

            for (int i = 0; i < numberOfMessages; i++)
            {
                int n = i;
                StringBuilder messageFromCore = new StringBuilder("                                                        ");
                MohidLandEngineDLLAccess.GetMessage(ref n, messageFromCore, (uint)messageFromCore.Length);
                message += "; ";
                message += messageFromCore.ToString().Trim();
            }
            throw new Exception(message);
        }

        /// <summary>
        /// Gets the name of the model
        /// </summary>
        /// <returns>Model Name</returns>
        public String GetModelID()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidLandEngineDLLAccess.GetModelID(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return stringBuilder.ToString().Trim();
        }

        /// <summary>
        /// Gets Start Instant of the Model
        /// </summary>
        public DateTime GetStartInstant()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidLandEngineDLLAccess.GetStartInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Stop Instant of the Model
        /// </summary>
        public DateTime GetStopInstant()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidLandEngineDLLAccess.GetStopInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Current Time of the Model
        /// </summary>
        public DateTime GetCurrentTime()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidLandEngineDLLAccess.GetCurrentInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Current Time Step of the Model
        /// </summary>
        public double GetCurrentTimeStep()
        {
            return MohidLandEngineDLLAccess.GetCurrentTimeStep();
        }

        /// <summary>
        /// Get number of Network Nodes
        /// </summary>
        public int GetNumberOfNodes()
        {
            return MohidLandEngineDLLAccess.GetNumberOfNodes();
        }



        public double GetXCoordinate(int nodeID)
        {
            return MohidLandEngineDLLAccess.GetXCoordinate(ref nodeID);
        }

        public double GetYCoordinate(int nodeID)
        {
            return MohidLandEngineDLLAccess.GetYCoordinate(ref nodeID);
        }

        private DateTime MohidTimeStringToDotNetTime(String mohidTimeString)
        {

            string[] TimeArray = mohidTimeString.Split(':');

            //Extracts year, month, day, hour, minute
            int year = Convert.ToInt32(TimeArray[0]);
            int month = Convert.ToInt32(TimeArray[1]);
            int day = Convert.ToInt32(TimeArray[2]);
            int hour = Convert.ToInt32(TimeArray[3]);
            int minute = Convert.ToInt32(TimeArray[4]);

            //Seconds
            int second = (int)System.Math.Floor(Convert.ToDouble(TimeArray[5]));

            //Milliseconds
            int millisecond = (int)(1000.0 * (Convert.ToDouble(TimeArray[5]) - (double)second));

            return new DateTime(year, month, day, hour, minute, second, millisecond);

            
        }


        public double GetFlowByNodeID(int nodeID)
        {
            return MohidLandEngineDLLAccess.GetFlowByNodeID(ref nodeID);
        }
    }
}
