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

        //TODO: Check how we can the correct instance IDs from the MOHID models. Implement a selector in the main, which receives the model ID?
        private int dischargeInstanceID = 1;
        private int horizontalGridInstanceID = 1;
        private int horizontalMapInstanceID = 1;
        private int hydrodynamicInstanceID = 1;


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

            SpatialReference spatialRef = new SpatialReference("ref");


            //Output exchange items - water level in each grid cell
            for (int i = 1; i <= MohidWaterEngine.GetIUB(horizontalGridInstanceID); i++)
            {
                for (int j = 1; j <= MohidWaterEngine.GetJUB(horizontalGridInstanceID); j++)
                {
                    if (MohidWaterEngine.IsWaterPoint(horizontalMapInstanceID, i, j))
                    {

                        String name = "i=" + i.ToString() + "/j=" + j.ToString();

                        ElementSet waterLevelPoint = new ElementSet("Grid Cell of MOHID Water", name, ElementType.XYPoint, spatialRef);
                        Element element = new Element(name);
                        element.AddVertex(new Vertex(MohidWaterEngine.GetCenterXCoordinate(horizontalGridInstanceID, i, j), MohidWaterEngine.GetCenterYCoordinate(horizontalGridInstanceID, i, j), 0));
                        waterLevelPoint.AddElement(element);

                        OutputExchangeItem waterLevel = new OutputExchangeItem();
                        waterLevel.Quantity = waterLevelQuantity;
                        waterLevel.ElementSet = waterLevelPoint;

                        outputExchangeItems.Add(waterLevel);

                        //DataOperation op = new DataOperation();
                        //op.AddArgument(new Argument("VerticalShift", "2.00", true, "Define numerical value (negative is towards earth centre)"));
                        //waterLevel.AddDataOperation(op);
                        

                    }
                }
            }


            //Flow input exchange to discharges configured as OpenMI Discharges
            for (int i = 1; i <= MohidWaterEngine.GetNumberOfDischarges(dischargeInstanceID); i++)
			{
                if (MohidWaterEngine.GetDischargeType(dischargeInstanceID, i) == 4)
                {

                    String pointName = MohidWaterEngine.GetDischargeName(dischargeInstanceID, i);

                    ElementSet dischargePoint = new ElementSet("Discharge Point of MOHID Water", i.ToString(), ElementType.XYPoint, new SpatialReference("ref"));

                    InputExchangeItem inputDischarge = new InputExchangeItem();
                    inputDischarge.Quantity = flowQuantity;

                    Element element = new Element("Point: " + pointName);
                    element.AddVertex(new Vertex(MohidWaterEngine.GetDischargeXCoordinate(dischargeInstanceID, i), MohidWaterEngine.GetDischargeYCoordinate(dischargeInstanceID, i), 0));
                    dischargePoint.AddElement(element);
                    
                    inputDischarge.ElementSet = dischargePoint;

                    inputExchangeItems.Add(inputDischarge);
                }
			}




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
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(MohidWaterEngine.GetCurrentTime().AddSeconds(MohidWaterEngine.GetCurrentTimeStep())));
        }

        public double GetMissingValueDefinition()
        {
            return -99;
        }

        public global::OpenMI.Standard.IValueSet GetValues(string QuantityID, string ElementSetID)
        {

            double[] returnValues;
            Char[] separator = new char[] { '/', '=' };

            if (QuantityID == "Water Level")
            {

                //Gets the water level at the specified point
                String[] substrings = ElementSetID.Split(separator);

                int i = Convert.ToInt32(substrings[1]);
                int j = Convert.ToInt32(substrings[3]);

                returnValues = new double[1];
                //returnValues[0] = MohidWaterEngine.GetWaterLevelAtPoint(hydrodynamicInstanceID, i, j);
                //TODO: Implement this vertical shift through data operation
                returnValues[0] = MohidWaterEngine.GetWaterLevelAtPoint(hydrodynamicInstanceID, i, j) - 2.00;
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
            if (QuantityID == "Flow")
            {
                double flow = ((ScalarSet)values).data[0];
                MohidWaterEngine.SetDischargeFlow(dischargeInstanceID, Convert.ToInt32(ElementSetID), flow);
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
