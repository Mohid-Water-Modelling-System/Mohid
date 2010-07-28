using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MOHID.OpenMI.Sdk.Wrapper;
using MOHID.OpenMI.Sdk.Backbone;
using MOHID.OpenMI.Sdk.DevelopmentSupport;
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
            Quantity inFlowQuantity = new Quantity(flowUnit, "description", "InFlow", global::OpenMI.Standard.ValueType.Scalar, flowDimension);

            int numberOfNodes = mohidLandEngine.GetNumberOfNodes();

            for (int i = 1; i <= numberOfNodes; i++)
            {
                OutputExchangeItem flowFromNode = new OutputExchangeItem();
                InputExchangeItem inFlowToNode = new InputExchangeItem();

                ElementSet node = new ElementSet("description", "Node: " + i.ToString(), ElementType.XYPoint, new SpatialReference("ref"));
                node.AddElement(new Element("Node:" + i.ToString()));
                node.Elements[0].AddVertex(new Vertex(mohidLandEngine.GetXCoordinate(i), mohidLandEngine.GetYCoordinate(i), 0));

                flowFromNode.ElementSet = node;
                flowFromNode.Quantity = flowQuantity;

                inFlowToNode.ElementSet = node;
                inFlowToNode.Quantity = inFlowQuantity;

                outputExchangeItems.Add(flowFromNode);
                inputExchangeItems.Add(inFlowToNode);
            }

        }


        public InputExchangeItem GetInputExchangeItem(int exchangeItemIndex)
        {
            throw new NotImplementedException();
        }

        public int GetInputExchangeItemCount()
        {
            throw new NotImplementedException();
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
            throw new NotImplementedException();
        }

        public int GetOutputExchangeItemCount()
        {
            throw new NotImplementedException();
        }

        public ITimeSpan GetTimeHorizon()
        {
            DateTime modelStart = mohidLandEngine.GetStartInstant();
            DateTime modelEnd = mohidLandEngine.GetStopInstant();

            double start = CalendarConverter.Gregorian2ModifiedJulian(modelStart);
            double end = CalendarConverter.Gregorian2ModifiedJulian(modelEnd);

            ITimeStamp tStart = new OpenMI.Sdk.Backbone.TimeStamp(start);
            ITimeStamp tEnd = new TimeStamp(end);

            return new MOHID.OpenMI.Sdk.Backbone.TimeSpan(tStart, tEnd);
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
            throw new NotImplementedException();
        }

        public string GetComponentID()
        {
            throw new NotImplementedException();
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
            throw new NotImplementedException();
        }

        public global::OpenMI.Standard.IValueSet GetValues(string QuantityID, string ElementSetID)
        {

            double[] returnValues;
            Char[] separator = new char[] { ':' };

            if (QuantityID == "Flow")
            {
                int index = Convert.ToInt32((ElementSetID.Split(separator))[1]);
                returnValues = new double[1];
                returnValues[0] = mohidLandEngine.GetFlowByNodeID(index);
            }
            else
            {
                throw new Exception("Illegal QuantityID in GetValues method in SimpleRiverEngine");
            }

            ScalarSet values = new ScalarSet(returnValues);
            return values;
        }


        public bool PerformTimeStep()
        {
            throw new NotImplementedException();
        }

        public void SetValues(string QuantityID, string ElementSetID, global::OpenMI.Standard.IValueSet values)
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
