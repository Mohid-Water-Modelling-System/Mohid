using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;

namespace MOHID.OpenMI.MohidWater.Wrapper
{

    public class MohidWaterEngineDotNetAccess
    {


        IntPtr _FortranDllHandle;

        /// <summary>
        /// Loads the Fortran Dll and initializes the model (calling the constructor)
        /// </summary>
        /// <param name="filePath">Path to the nomfich.dat file</param>
        public void Initialize(string filePath)
        {
            //Loads the library
            _FortranDllHandle =  Kernel32Wrapper.LoadLibrary(@"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidWaterEngine\Release OpenMI\MohidWaterEngine.dll");

            //Sets the directory temporary to the exe dir of the model
            String currentDir = Environment.CurrentDirectory;
            Environment.CurrentDirectory = System.IO.Path.GetDirectoryName(filePath);

            //Calls the constructor and reads data files, etc
            if (!(MohidWaterEngineDLLAccess.Initialize(filePath, ((uint)filePath.Length))))
            {
                CreateAndThrowException();
            }

            Environment.CurrentDirectory = currentDir;
        
        }

        /// <summary>
        /// Performs a single time step
        /// </summary>
        public void PerformTimeStep()
        {
            if (!(MohidWaterEngineDLLAccess.PerformTimeStep()))
            {
                CreateAndThrowException();
            }
        }

        /// <summary>
        /// Calls the model destructor (closes data files)
        /// </summary>
        public void Finish()
        {
            if (!(MohidWaterEngineDLLAccess.Finish()))
            {
                CreateAndThrowException();
            }
            while (Kernel32Wrapper.FreeLibrary(_FortranDllHandle)) ;
        }

        public void Dispose()
        {
            if (!MohidWaterEngineDLLAccess.Dispose())
            {
                CreateAndThrowException();
            }
        }

        /// <summary>
        /// Runs the whole working cycle once - Testing Only
        /// </summary>
        public void RunSimulation()
        {
            if (!(MohidWaterEngineDLLAccess.RunSimulation()))
            {
                CreateAndThrowException();
            }
        }

        private void CreateAndThrowException()
        {
            int numberOfMessages = 0;
            numberOfMessages = MohidWaterEngineDLLAccess.GetNumberOfMessages();
            string message = "Error Messages from MOHID Water Engine";

            for (int i = 0; i < numberOfMessages; i++)
            {
                int n = i;
                StringBuilder messageFromCore = new StringBuilder("                                                        ");
                MohidWaterEngineDLLAccess.GetMessage(ref n, messageFromCore, (uint)messageFromCore.Length);
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
            if (!MohidWaterEngineDLLAccess.GetModelID(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return stringBuilder.ToString().Trim();
        }

        /// <summary>
        /// Gets Start Instant of the Model
        /// </summary>
        public DateTime GetStartInstant()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidWaterEngineDLLAccess.GetStartInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Stop Instant of the Model
        /// </summary>
        public DateTime GetStopInstant()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidWaterEngineDLLAccess.GetStopInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Current Time of the Model
        /// </summary>
        public DateTime GetCurrentTime()
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidWaterEngineDLLAccess.GetCurrentInstant(stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return MohidTimeStringToDotNetTime(stringBuilder.ToString().Trim());
        }

        /// <summary>
        /// Gets Current Time Step of the Model
        /// </summary>
        public double GetCurrentTimeStep()
        {
            return MohidWaterEngineDLLAccess.GetCurrentTimeStep();
        }

        #region Module Discharges

        ///// <summary>
        ///// Get number of discharges
        ///// </summary>
        public int GetNumberOfDischarges(int dischargeInstanceID)
        {
            return MohidWaterEngineDLLAccess.GetNumberOfDischarges(ref dischargeInstanceID);
        }

        public int GetDischargeType(int dischargeInstanceID, int dischargeID)
        {
            return MohidWaterEngineDLLAccess.GetDischargeType(ref dischargeInstanceID, ref dischargeID);
        }

        public double GetDischargeXCoordinate(int dischargeInstanceID, int dischargeID)
        {
            return MohidWaterEngineDLLAccess.GetDischargeXCoordinate(ref dischargeInstanceID, ref dischargeID);
        }

        public double GetDischargeYCoordinate(int dischargeInstanceID, int dischargeID)
        {
            return MohidWaterEngineDLLAccess.GetDischargeYCoordinate(ref dischargeInstanceID, ref dischargeID);
        }

        public string GetDischargeName(int dischargeInstanceID, int dischargeID)
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidWaterEngineDLLAccess.GetDischargeName(ref dischargeInstanceID, ref dischargeID, stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return stringBuilder.ToString().Trim();
        }

        public void SetDischargeFlow(int dischargeInstanceID, int dischargeID, double flow)
        {
            MohidWaterEngineDLLAccess.SetDischargeFlow(ref dischargeInstanceID, ref dischargeID, ref flow);
        }

        public int GetNumberOfDischargeProperties(int dischargeInstanceID, int dischargeID)
        {
            return MohidWaterEngineDLLAccess.GetNumberOfDischargeProperties(ref dischargeInstanceID, ref dischargeID);
        }

        public int GetDischargePropertyID(int dischargeInstanceID, int dischargeID, int idx)
        {
            return MohidWaterEngineDLLAccess.GetDischargePropertyID(ref dischargeInstanceID, ref dischargeID, ref idx);
        }

        public void SetDischargeConcentration(int dischargeInstanceID, int dischargeID, int propertyID, double concentration)
        {
            MohidWaterEngineDLLAccess.SetDischargeConcentration(ref dischargeInstanceID, ref dischargeID, ref propertyID, ref concentration);
        }


        #endregion

        #region Module HorinzontalGrid

        public int GetIUB(int horizontalGridInstanceID)
        {
            return MohidWaterEngineDLLAccess.GetIUB(ref horizontalGridInstanceID);
        }

        public int GetJUB(int horizontalGridInstanceID)
        {
            return MohidWaterEngineDLLAccess.GetJUB(ref horizontalGridInstanceID);
        }

        public bool IsWaterPoint(int horizontalGridInstanceID, int i, int j)
        {
            return MohidWaterEngineDLLAccess.IsWaterPoint(ref horizontalGridInstanceID, ref i, ref j);
        }

        public double GetCenterXCoordinate(int horizontalGridInstanceID, int i, int j)
        {
            return MohidWaterEngineDLLAccess.GetCenterXCoordinate(ref horizontalGridInstanceID, ref i, ref j);
        }

        public double GetCenterYCoordinate(int horizontalGridInstanceID, int i, int j)
        {
            return MohidWaterEngineDLLAccess.GetCenterYCoordinate(ref horizontalGridInstanceID, ref i, ref j);
        }

        public void GetGridCellCoordinates(int horizontalGridInstanceID, int i, int j, ref double[] xCoords, ref double[] yCoords)
        {
            MohidWaterEngineDLLAccess.GetGridCellCoordinates(ref horizontalGridInstanceID, ref i, ref j, xCoords, yCoords);
        }

        #endregion

        #region Module Hydrodynamic

        public double GetWaterLevelAtPoint(int hydrodynamicInstanceID, int i, int j)
        {
            return MohidWaterEngineDLLAccess.GetWaterLevelAtPoint(ref hydrodynamicInstanceID, ref i, ref j);
        }

        public void GetWaterLevel1D(int hydrodynamicInstanceID, int numberOfComputePoints, ref double[] waterLevels1D)
        {
            MohidWaterEngineDLLAccess.GetWaterLevel1D(ref hydrodynamicInstanceID, ref numberOfComputePoints, waterLevels1D);
        }

        #endregion

        #region Module Waterproperties

        public int GetNumberOfProperties(int waterPropertiesInstanceID)
        {
            return MohidWaterEngineDLLAccess.GetNumberOfProperties(ref waterPropertiesInstanceID);
        }

        public int GetPropertyIDNumber(int waterPropertiesInstanceID, int idx)
        {
            return MohidWaterEngineDLLAccess.GetWaterPropertiesPropertyID(ref waterPropertiesInstanceID, ref idx);
        }

        public string GetPropertyNameByIDNumber(int propertyID)
        {
            StringBuilder stringBuilder = new StringBuilder("                         ");
            if (!MohidWaterEngineDLLAccess.GetPropertyNameByID(ref propertyID, stringBuilder, (uint)stringBuilder.Length))
                CreateAndThrowException();
            return stringBuilder.ToString().Trim();
        }

        public double GetConcentrationAtPoint(int waterPropertiesInstanceID, int propertyID, int i, int j)
        {
            return MohidWaterEngineDLLAccess.GetConcentrationAtPoint(ref waterPropertiesInstanceID, ref propertyID, ref i, ref j);
        }

        public void GetConcentration1D(int waterPropertiesInstanceID, int propertyID, int numberOfcomputePoints, ref double[] concentration1D)
        {
            MohidWaterEngineDLLAccess.GetConcentration1D(ref waterPropertiesInstanceID, ref propertyID,
                                                                ref numberOfcomputePoints, concentration1D);
        }

        #endregion

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



    }
}
