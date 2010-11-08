using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace MOHID.OpenMI.MohidLand.Wrapper
{
    public class MohidLandEngineDLLAccess
    {
        //TODO: Check how to set this path during runtime ou by compiler reference...
        private const string dllPath = @"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidLandEngine\Release OpenMI\MohidLandEngine.dll";

        /// <summary>
        /// Calls the Initialize public method in the MOHID Land DLL
        /// </summary>
        /// <param name="filePath">Path to the nomfich file</param>
        /// <param name="length">Length of the path</param>
        /// <returns>TRUE if successful, FALSE otherwise</returns>
        [DllImport(dllPath, 
            EntryPoint = "INITIALIZE", 
            SetLastError = true, 
            ExactSpelling = true, 
            CallingConvention = CallingConvention.Cdecl)]
        public static extern bool Initialize(string filePath, uint length);

        [DllImport(dllPath, EntryPoint = "PERFORMTIMESTEP", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool PerformTimeStep();

        [DllImport(dllPath, EntryPoint = "GETSTARTINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetStartInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETSTOPINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetStopInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETCURRENTINSTANT", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetCurrentInstant([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        /// <summary>
        /// Gets the current time step of the MOHID Land model
        /// </summary>
        /// <returns></returns>
        [DllImport(dllPath, EntryPoint = "GETCURRENTTIMESTEP", 
            SetLastError = true, 
            ExactSpelling = true, 
            CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetCurrentTimeStep();


        [DllImport(dllPath, EntryPoint = "GETMODELID", SetLastError = true, ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetModelID([MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFMESSAGES",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfMessages();

        [DllImport(dllPath,EntryPoint = "GETMESSAGE",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetMessage(ref int messageID, [MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        [DllImport(dllPath,EntryPoint = "GETPROPERTYNAMEBYID",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetPropertyNameByID(ref int propertyID, [MarshalAs(UnmanagedType.LPStr)] StringBuilder id, uint length);

        #region Module Drainage Network

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFNODES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfNodes(ref int drainageNetworkInstanceID);

        [DllImport(dllPath, EntryPoint = "GETXCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetXCoordinate(ref int drainageNetworkInstanceID, ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETYCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetYCoordinate(ref int drainageNetworkInstanceID, ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETFLOWBYNODEID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetFlowByNodeID(ref int drainageNetworkInstanceID, ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETOUTLETFLOW", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetOutletFlow(ref int drainageNetworkInstanceID);

        [DllImport(dllPath, EntryPoint = "SETDOWNSTREAMWATERLEVEL", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool SetDownStreamWaterLevel(ref int drainageNetworkInstanceID, ref double waterLevel);

        [DllImport(dllPath, EntryPoint = "GETOUTLETNODEID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetOutletNodeID(ref int drainageNetworkInstanceID);

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFPROPERTIES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfProperties(ref int drainageNetworkInstanceID);

        [DllImport(dllPath, EntryPoint = "GETDRAINAGENETWORKPROPERTYID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetDrainageNetworkPropertyID(ref int drainageNetworkInstanceID, ref int idx);

        [DllImport(dllPath, EntryPoint = "GETOUTLETFLOWCONCENTRATION", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetOutletFlowConcentration(ref int drainageNetworkInstanceID, ref int propertyID);

        [DllImport(dllPath, EntryPoint = "SETDOWNSTREAMCONCENTRATION", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool SetDownStreamConcentration(ref int drainageNetworkInstanceID, ref int propertyID, ref double concentration);


        #endregion

        [DllImport(dllPath,EntryPoint = "RUNSIMULATION",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool RunSimulation();

        #region Destructor

        [DllImport(dllPath,EntryPoint = "FINISH",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        public static extern bool Finish();

        //[DllImport(dllPath,EntryPoint = "DISPOSE",SetLastError = true,ExactSpelling = true,CallingConvention = CallingConvention.Cdecl)]
        //public static extern bool Dispose();

        
        #endregion

    }
}


