using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace MOHID.OpenMI.MohidWater.Wrapper
{
    public class MohidWaterEngineDLLAccess
    {
        //TODO: Check how to set this path during runtime ou by compiler reference...
        //private const string dllPath = @"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidWaterEngine\Release OpenMI\MohidWaterEngine.dll";
        private const string dllPath = @"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidWaterEngine\Debug OpenMI\MohidWaterEngine.dll";

        [DllImport(dllPath, EntryPoint = "INITIALIZE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool Initialize(string filePath, uint length);

        [DllImport(dllPath, EntryPoint = "PERFORMTIMESTEP", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool PerformTimeStep();

        [DllImport(dllPath, EntryPoint = "GETSTARTINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetStartInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETSTOPINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetStopInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETCURRENTINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetCurrentInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETCURRENTTIMESTEP", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetCurrentTimeStep();

        [DllImport(dllPath, EntryPoint = "GETMODELID", SetLastError = true, ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetModelID([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFMESSAGES",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfMessages();

        [DllImport(dllPath,EntryPoint = "GETMESSAGE",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetMessage(ref int messageID, [MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETPROPERTYNAMEBYID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetPropertyNameByID(ref int propertyID, [MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        #region Module Discharges

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFDISCHARGES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfDischarges(ref int instanceID);

        [DllImport(dllPath, EntryPoint = "GETDISCHARGETYPE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetDischargeType(ref int instanceID, ref int dischargeID);

        [DllImport(dllPath, EntryPoint = "GETDISCHARGEXCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetDischargeXCoordinate(ref int instanceID, ref int dischargeID);

        [DllImport(dllPath, EntryPoint = "GETDISCHARGEYCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetDischargeYCoordinate(ref int instanceID, ref int dischargeID);

        [DllImport(dllPath, EntryPoint = "GETDISCHARGENAME", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetDischargeName(ref int instanceID, ref int dischargeID, [MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "SETDISCHARGEFLOW", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool SetDischargeFlow(ref int instanceID, ref int dischargeID, ref double flow);

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFDISCHARGEPROPERTIES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfDischargeProperties(ref int instanceID, ref int dischargeID);

        [DllImport(dllPath, EntryPoint = "GETDISCHARGEPROPERTYID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetDischargePropertyID(ref int instanceID, ref int dischargeID, ref int idx);

        [DllImport(dllPath, EntryPoint = "SETDISCHARGECONCENTRATION", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool SetDischargeConcentration(ref int instanceID, ref int dischargeID, ref int propertyID, ref double concentration);

        #endregion

        #region Module HorizontalGrid / Map

        [DllImport(dllPath, EntryPoint = "GETIUB", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetIUB(ref int horizontalGridInstanceID);

        [DllImport(dllPath, EntryPoint = "GETJUB", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetJUB(ref int horizontalGridInstanceID);

        [DllImport(dllPath, EntryPoint = "ISWATERPOINT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool IsWaterPoint(ref int horizontalGridInstanceID, ref int i, ref int j);

        [DllImport(dllPath, EntryPoint = "GETCENTERXCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetCenterXCoordinate(ref int horizontalGridInstanceID, ref int i, ref int j);

        [DllImport(dllPath, EntryPoint = "GETCENTERYCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetCenterYCoordinate(ref int horizontalGridInstanceID, ref int i, ref int j);

        [DllImport(dllPath, EntryPoint = "GETGRIDCELLCOORDINATES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetGridCellCoordinates(ref int horizontalGridInstanceID, ref int i, ref int j, double[] xCoords, double[] yCoords);

        #endregion

        #region

        [DllImport(dllPath, EntryPoint = "GETWATERLEVELATPOINT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetWaterLevelAtPoint(ref int hydrodynamicInstanceID, ref int i, ref int j);

        [DllImport(dllPath, EntryPoint = "GETWATERLEVEL1D", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetWaterLevel1D(ref int hydrodynamicInstanceID, ref int numberOfWaterPoints, double[] waterlevels1D);
        

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFPROPERTIES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfProperties(ref int waterPropertiesID);

        [DllImport(dllPath, EntryPoint = "GETWATERPROPERTIESPROPERTYID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetWaterPropertiesPropertyID(ref int waterPropertiesID, ref int idx);

        [DllImport(dllPath, EntryPoint = "GETCONCENTRATIONATPOINT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetConcentrationAtPoint(ref int waterPropertiesID, ref int propertyID, ref int i, ref int j);

        [DllImport(dllPath, EntryPoint = "GETCONCENTRATION1D", SetLastError = true, ExactSpelling = true,
            CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetConcentration1D(ref int waterPropertiesID, ref int propertyID,
                                                       ref int numberOfWaterPoints, double[] concentration1D);

        #endregion

        [DllImport(dllPath,EntryPoint = "RUNSIMULATION",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool RunSimulation();

        #region Destructor

        [DllImport(dllPath,EntryPoint = "FINISH",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool Finish();

        [DllImport(dllPath,EntryPoint = "DISPOSE",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool Dispose();

        
        #endregion




    }
}


