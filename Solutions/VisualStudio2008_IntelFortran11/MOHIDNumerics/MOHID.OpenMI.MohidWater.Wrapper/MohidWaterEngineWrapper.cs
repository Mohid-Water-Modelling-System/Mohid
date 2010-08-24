using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Oatc.OpenMI.Sdk.Wrapper;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.DevelopmentSupport;
using OpenMI.Standard;
using System.Collections;

namespace MOHID.OpenMI.MohidWater.Wrapper
{
    public class MohidWaterEngineWrapper : IEngine
    {

        #region Fields

        private MohidWaterEngineDotNetAccess MohidWaterEngine;
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
            MohidWaterEngine = new MohidWaterEngineDotNetAccess();
            MohidWaterEngine.Initialize(properties["FilePath"].ToString());

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
            //outletNode.AddElement(new Element("Outlet"));
            //int outletNodeID = MohidWaterEngine.GetOutletNodeID();
            //outletNode.Elements[0].AddVertex(new Vertex(MohidWaterEngine.GetXCoordinate(outletNodeID), MohidWaterEngine.GetYCoordinate(outletNodeID), 0));
            //outletFlow.ElementSet = outletNode;
            //outletFlow.Quantity = flowQuantity;

            outputExchangeItems.Add(outletFlow);

            //Discharges
            InputExchangeItem inputDischarge = new InputExchangeItem();
            inputDischarge.Quantity = flowQuantity;
            inputDischarge.ElementSet = outletNode;

            inputExchangeItems.Add(inputDischarge);


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
            return MohidWaterEngine.GetModelID();
        }

        public string GetModelID()
        {
            return MohidWaterEngine.GetModelID();
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
            DateTime modelStart = MohidWaterEngine.GetStartInstant();
            DateTime modelEnd = MohidWaterEngine.GetStopInstant();

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
            MohidWaterEngine.Dispose();
            MohidWaterEngine = null;
        }

        public void Finish()
        {
            MohidWaterEngine.Finish();
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
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(MohidWaterEngine.GetCurrentTime()));
        }

        public global::OpenMI.Standard.ITimeStamp GetEarliestNeededTime()
        {
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(MohidWaterEngine.GetCurrentTime().AddSeconds(MohidWaterEngine.GetCurrentTimeStep())));
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
                //returnValues[0] = MohidWaterEngine.GetOutletFlow();
            }
            else
            {
                throw new Exception("Illegal QuantityID in GetValues method in MohidWaterEngineWrapper");
            }

            ScalarSet values = new ScalarSet(returnValues);
            return values;
        }

        public void SetValues(string QuantityID, string ElementSetID, global::OpenMI.Standard.IValueSet values)
        {
            if (QuantityID == "Water Level")
            {
                double waterLevel = ((ScalarSet)values).data[0];
                //MohidWaterEngine.SetDownstreamWaterLevel(waterLevel);
            }
            else
            {
                throw new Exception("Illegal QuantityID in SetValues method in MohidWaterEngineWrapper");
            }
        }


        public bool PerformTimeStep()
        {
            MohidWaterEngine.PerformTimeStep();
            return true;
        }

 
        #endregion
    }
}
