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
        private const string dllPath = @"D:\Software\Mohid\MOHID.Numerics\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics\MohidLandEngine\Debug\MohidLandEngine.dll";

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

        [DllImport(dllPath, EntryPoint = "GETNUMBEROFNODES", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfNodes();

        [DllImport(dllPath, EntryPoint = "GETXCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetXCoordinate(ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETYCOORDINATE", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetYCoordinate(ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETFLOWBYNODEID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetFlowByNodeID(ref int nodeID);

        [DllImport(dllPath, EntryPoint = "GETOUTLETFLOW", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GetOutletFlow();

        [DllImport(dllPath, EntryPoint = "SETDOWNSTREAMWATERLEVEL", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern bool SetDownStreamWaterLevel();

        [DllImport(dllPath, EntryPoint = "GETOUTLETNODEID", SetLastError = true, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetOutletNodeID();

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


