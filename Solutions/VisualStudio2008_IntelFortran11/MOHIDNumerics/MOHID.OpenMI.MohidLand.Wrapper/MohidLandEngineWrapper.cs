using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Oatc.OpenMI.Sdk.Wrapper;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.DevelopmentSupport;
using OpenMI.Standard;
using System.Collections;

namespace MOHID.OpenMI.MohidLand.Wrapper
{
    public class MohidLandEngineWrapper : IEngine
    {

        #region Fields

        private MohidLandEngineDotNetAccess mohidLandEngine;
        private ArrayList inputExchangeItems;
        private ArrayList outputExchangeItems;

        #endregion

        #region IEngine Members

        public void Initialize(System.Collections.Hashtable properties)
        {
            //List of exchange Items
            inputExchangeItems = new ArrayList();
            outputExchangeItems = new ArrayList();

            //Initializes Engine
            mohidLandEngine = new MohidLandEngineDotNetAccess();
            mohidLandEngine.Initialize(properties["FilePath"].ToString());

            // -- Build exchange items ---
            Dimension flowDimension = new Dimension();
            Unit flowUnit = new Unit("m3/sec", 1, 0, "m3/sec");
            Quantity flowQuantity = new Quantity(flowUnit, "description", "Flow", global::OpenMI.Standard.ValueType.Scalar, flowDimension);

            Dimension waterlevelDimension = new Dimension();
            Unit waterlevelUnit = new Unit("m", 1, 0, "m");
            Quantity waterLevelQuantity = new Quantity(waterlevelUnit, "description", "Water Level", global::OpenMI.Standard.ValueType.Scalar, waterlevelDimension);


            //Flow at the outlet
            OutputExchangeItem outletFlow = new OutputExchangeItem();
            outletFlow.Quantity = flowQuantity;
            ElementSet outletNode = new ElementSet("description", "Outlet", ElementType.XYPoint, new SpatialReference("ref"));
            outletNode.AddElement(new Element("Outlet"));
            int outletNodeID = mohidLandEngine.GetOutletNodeID();
            outletNode.Elements[0].AddVertex(new Vertex(mohidLandEngine.GetXCoordinate(outletNodeID), mohidLandEngine.GetYCoordinate(outletNodeID), 0));
            outletFlow.ElementSet = outletNode;
            outletFlow.Quantity = flowQuantity;

            outputExchangeItems.Add(outletFlow);

            //Water level at the outlet
            InputExchangeItem outletLevel = new InputExchangeItem();
            outletLevel.Quantity = waterLevelQuantity;
            outletLevel.ElementSet = outletNode;
            outletLevel.Quantity = waterLevelQuantity;

            inputExchangeItems.Add(outletLevel);


        }


        public InputExchangeItem GetInputExchangeItem(int exchangeItemIndex)
        {
            return this.inputExchangeItems[exchangeItemIndex] as InputExchangeItem;
        }

        public int GetInputExchangeItemCount()
        {
            return this.inputExchangeItems.Count;
        }

        public string GetModelDescription()
        {
            return mohidLandEngine.GetModelID();
        }

        public string GetModelID()
        {
            return mohidLandEngine.GetModelID();
        }

        public OutputExchangeItem GetOutputExchangeItem(int exchangeItemIndex)
        {
            return this.outputExchangeItems[exchangeItemIndex] as OutputExchangeItem;
        }

        public int GetOutputExchangeItemCount()
        {
            return this.outputExchangeItems.Count;
        }

        public ITimeSpan GetTimeHorizon()
        {
            DateTime modelStart = mohidLandEngine.GetStartInstant();
            DateTime modelEnd = mohidLandEngine.GetStopInstant();

            double start = CalendarConverter.Gregorian2ModifiedJulian(modelStart);
            double end = CalendarConverter.Gregorian2ModifiedJulian(modelEnd);

            ITimeStamp tStart = new Oatc.OpenMI.Sdk.Backbone.TimeStamp(start);
            ITimeStamp tEnd = new TimeStamp(end);

            return new Oatc.OpenMI.Sdk.Backbone.TimeSpan(tStart, tEnd);
        }

        #endregion

        #region IRunEngine Members

        public void Dispose()
        {
            mohidLandEngine.Dispose();
            mohidLandEngine = null;
        }

        public void Finish()
        {
            mohidLandEngine.Finish();
        }

        public string GetComponentDescription()
        {
            return "GetComponentDescription";
        }

        public string GetComponentID()
        {
            return "GetComponentID";
        }

        public global::OpenMI.Standard.ITime GetCurrentTime()
        {
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(mohidLandEngine.GetCurrentTime()));
        }

        public global::OpenMI.Standard.ITimeStamp GetEarliestNeededTime()
        {
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(mohidLandEngine.GetCurrentTime().AddSeconds(mohidLandEngine.GetCurrentTimeStep())));
        }

        public global::OpenMI.Standard.ITime GetInputTime(string QuantityID, string ElementSetID)
        {
            throw new NotImplementedException();
        }

        public double GetMissingValueDefinition()
        {
            return -99;
        }

        public global::OpenMI.Standard.IValueSet GetValues(string QuantityID, string ElementSetID)
        {

            double[] returnValues;
            Char[] separator = new char[] { ':' };

            if (QuantityID == "Flow")
            {
                returnValues = new double[1];
                returnValues[0] = mohidLandEngine.GetOutletFlow();
            }
            else
            {
                throw new Exception("Illegal QuantityID in GetValues method in MohidLandEngineWrapper");
            }

            ScalarSet values = new ScalarSet(returnValues);
            return values;
        }

        public void SetValues(string QuantityID, string ElementSetID, global::OpenMI.Standard.IValueSet values)
        {
            if (QuantityID == "Water Level")
            {
                double waterLevel = ((ScalarSet)values).data[0];
                mohidLandEngine.SetDownstreamWaterLevel(waterLevel);
            }
            else
            {
                throw new Exception("Illegal QuantityID in SetValues method in MohidLandEngineWrapper");
            }
        }


        public bool PerformTimeStep()
        {
            mohidLandEngine.PerformTimeStep();
            return true;
        }

 
        #endregion
    }
}
