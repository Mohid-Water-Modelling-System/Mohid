using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Oatc.OpenMI.Sdk.Wrapper;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.DevelopmentSupport;
using OpenMI.Standard;
using System.Collections;
using ValueType = OpenMI.Standard.ValueType;
using Oatc.OpenMI.Sdk.Spatial;
using System.Diagnostics;

namespace MOHID.OpenMI.MohidLand.Wrapper
{
    public class MohidLandEngineWrapper : IEngine
    {

        #region Fields

        private MohidLandEngineDotNetAccess mohidLandEngine;
        private IList<InputExchangeItem> inputExchangeItems;
        private IList<OutputExchangeItem> outputExchangeItems;

        //TODO: Check how we can the correct instance IDs from the MOHID models. Implement a selector in the main, which receives the model ID?
        private int drainageNetworkInstanceID = 1;
        private int horizontalGridInstanceID = 1;
        private int horizontalMapInstanceID = 1;
        private int runoffInstanceID = 1;

        private Quantity qtdOutflow;            //Flow at the outlet                                     -> output
        private Quantity qtdOutletLevel;        //Level at the outlet                                    -> input
        private Quantity qtdWaterColumn;        //Overland Water Column                                  -> output
        private Quantity qtdDischarges;         //Discharges into the sewer system from the water column -> input
        private Quantity qtdFlowToStorm;        //Flow from nodes to sewer system                        -> output
        private Quantity qtdFlowFromStrom;      //Flow from sewer system to nodes                        -> input

        private int numberOfWaterPoints;

        private Stopwatch getValuesWatch;
        private Stopwatch setValuesWatch;
        private Stopwatch performStepWatch;
 
        #endregion

        #region IEngine Members

        public void Initialize(System.Collections.Hashtable properties)
        {

            getValuesWatch = new Stopwatch();
            setValuesWatch = new Stopwatch();
            performStepWatch = new Stopwatch();

            //List of exchange Items
            inputExchangeItems = new List<InputExchangeItem>();
            outputExchangeItems = new List<OutputExchangeItem>();

            //Initializes Engine
            mohidLandEngine = new MohidLandEngineDotNetAccess();
            mohidLandEngine.Initialize(properties["FilePath"].ToString());


            //
            //Dimensions
            //
            Dimension dimFlow = new Dimension();
            dimFlow.SetPower(DimensionBase.Length, 3);
            dimFlow.SetPower(DimensionBase.Time, -1);

            Dimension dimWaterlevel = new Dimension();
            dimWaterlevel.SetPower(DimensionBase.Length, 1);

            Dimension dimWaterColumn = new Dimension();
            dimWaterColumn.SetPower(DimensionBase.Length, 1);

            Dimension dimConcentration = new Dimension();
            dimConcentration.SetPower(DimensionBase.Mass, 1);
            dimConcentration.SetPower(DimensionBase.Length, -3);

            //
            //Units
            //
            Unit unitFlow = new Unit("m3/sec", 1, 0, "cubic meter per second");
            Unit unitWaterLevel = new Unit("m", 1, 0, "meters above mean sea level");
            Unit unitWaterColumn = new Unit("m", 1, 0, "meters above ground");
            Unit unitConcentration = new Unit("mg/l", 1, 0, "miligram per liter");

            //
            //Quantities
            //
            qtdOutflow = new Quantity(unitFlow, "Flow discharge at the outlet", "Outlet Flow", ValueType.Scalar, dimFlow);
            qtdOutletLevel = new Quantity(unitWaterLevel, "Waterlevel at the outlet", "OutletLevel",
                                                   ValueType.Scalar, dimWaterlevel);
            qtdWaterColumn = new Quantity(unitWaterColumn, "Ponded Water Column", "WaterColumn",
                                                   ValueType.Scalar, dimWaterColumn);
            qtdDischarges = new Quantity(unitFlow, "Distributed discharges (sewer sinks)", "Discharges", ValueType.Scalar,
                                                  dimFlow);

            qtdFlowToStorm = new Quantity(unitFlow, "Flow from the network to the storm water system (inlets)",
                                          "Storm Water Out Flow", ValueType.Scalar, dimFlow);

            qtdFlowFromStrom = new Quantity(unitFlow, "Flow from the storm water system to the network (discharges)",
                              "Storm Water In Flow", ValueType.Scalar, dimFlow);

            //
            //Spatial Reference
            //
            SpatialReference spatialReference = new SpatialReference("spatial reference");

            //
            //Element Sets
            //

            //Model Grid
            ElementSet modelGrid = new ElementSet("Model Grid Points of all Compute Points", "Model Grid",
                                       ElementType.XYPolygon, spatialReference);

            //Output exchange items - properties in each grid cell (surface only)
            numberOfWaterPoints = 0;
            for (int j = 1; j <= mohidLandEngine.GetJUB(horizontalGridInstanceID); j++)
            {
                for (int i = 1; i <= mohidLandEngine.GetIUB(horizontalGridInstanceID); i++)
                {
                    if (mohidLandEngine.IsWaterPoint(horizontalMapInstanceID, i, j))
                    {
                        String name = "i=" + i.ToString() + "/j=" + j.ToString();

                        double[] xCoords = new double[5];
                        double[] yCoords = new double[5];
                        mohidLandEngine.GetGridCellCoordinates(horizontalGridInstanceID, i, j, ref xCoords, ref yCoords);

                        Element element = new Element(name);
                        element.AddVertex(new Vertex(xCoords[0], yCoords[0], 0));
                        element.AddVertex(new Vertex(xCoords[1], yCoords[1], 0));
                        element.AddVertex(new Vertex(xCoords[2], yCoords[2], 0));
                        element.AddVertex(new Vertex(xCoords[3], yCoords[3], 0));

                        modelGrid.AddElement(element);

                        numberOfWaterPoints++;

                    }
                }
            }

            //Outlet Node
            ElementSet outletNode = new ElementSet("Outlet node", "Outlet", ElementType.XYPoint, spatialReference);
            Element outletElement = new Element("Outlet");
            int outletNodeID = mohidLandEngine.GetOutletNodeID(drainageNetworkInstanceID);
            outletElement.AddVertex(new Vertex(mohidLandEngine.GetXCoordinate(drainageNetworkInstanceID, outletNodeID),
                                               mohidLandEngine.GetYCoordinate(drainageNetworkInstanceID, outletNodeID),
                                               0));
            outletNode.AddElement(outletElement);

            //Outflow to Storm Water Model
            ElementSet stormWaterOutflowNodes = new ElementSet("Nodes which provide flow to the Storm Water System",
                                                          "Storm Water Inlets", ElementType.XYPoint, spatialReference);
            int numberOfOutflowNodes = mohidLandEngine.GetNumberOfStormWaterOutFlowNodes(drainageNetworkInstanceID);
            int[] outflowNodeIDs = new int[numberOfOutflowNodes];
            mohidLandEngine.GetStormWaterOutflowIDs(drainageNetworkInstanceID, numberOfOutflowNodes, ref outflowNodeIDs);
            for (int i = 1; i <= numberOfOutflowNodes; i++)
            {
                int nodeID = outflowNodeIDs[i - 1];

                Element element = new Element(nodeID.ToString());
                element.AddVertex(new Vertex(mohidLandEngine.GetXCoordinate(drainageNetworkInstanceID, nodeID),
                                             mohidLandEngine.GetYCoordinate(drainageNetworkInstanceID, nodeID),
                                             0));
                stormWaterOutflowNodes.AddElement(element);
            }

            //Inflow from Storm Water Model
            ElementSet stormWaterInflowNodes = new ElementSet("Nodes which receive flow to the Storm Water System",
                                                          "Storm Water Outlets", ElementType.XYPoint, spatialReference);
            int numberOfInflowNodes = mohidLandEngine.GetNumberOfStormWaterInFlowNodes(drainageNetworkInstanceID);
            if (numberOfInflowNodes > 0)
            {
                int[] inflowNodeIDs = new int[numberOfInflowNodes];
                mohidLandEngine.GetStormWaterInflowIDs(drainageNetworkInstanceID, numberOfOutflowNodes,
                                                       ref inflowNodeIDs);
                for (int i = 1; i <= numberOfInflowNodes; i++)
                {
                    int nodeID = inflowNodeIDs[i - 1];

                    Element element = new Element(nodeID.ToString());
                    element.AddVertex(new Vertex(mohidLandEngine.GetXCoordinate(drainageNetworkInstanceID, nodeID),
                                                 mohidLandEngine.GetYCoordinate(drainageNetworkInstanceID, nodeID),
                                                 0));
                    stormWaterInflowNodes.AddElement(element);
                }
            }


            //
            //Output Exchange Items
            //

            //Flow at the outlet
            OutputExchangeItem outletFlow = new OutputExchangeItem();
            outletFlow.Quantity = qtdOutflow;
            outletFlow.ElementSet = outletNode;
            outputExchangeItems.Add(outletFlow);

            //Overland water column
            OutputExchangeItem overlandWaterColumn = new OutputExchangeItem();
            overlandWaterColumn.Quantity = qtdWaterColumn;
            overlandWaterColumn.ElementSet = modelGrid;
            outputExchangeItems.Add(overlandWaterColumn);

            //Flow to the Storm Water Model
            if (stormWaterOutflowNodes.ElementCount > 0)
            {
                OutputExchangeItem stormWaterOutFlow = new OutputExchangeItem();
                stormWaterOutFlow.Quantity = qtdFlowToStorm;
                stormWaterOutFlow.ElementSet = stormWaterOutflowNodes;
                outputExchangeItems.Add(stormWaterOutFlow);
            }

            //
            //Input Exchange Items
            //

            //Water level at the outlet
            InputExchangeItem outletLevel = new InputExchangeItem();
            outletLevel.Quantity = qtdOutletLevel;
            outletLevel.ElementSet = outletNode;
            inputExchangeItems.Add(outletLevel);

            //Distributed discharges
            InputExchangeItem dischargeInflow = new InputExchangeItem();
            dischargeInflow.Quantity = qtdDischarges;
            dischargeInflow.ElementSet = modelGrid;
            inputExchangeItems.Add(dischargeInflow);

            //Flow from the Storm Water Model
            if (stormWaterInflowNodes.ElementCount > 0)
            {
                InputExchangeItem stormWaterInFlow = new InputExchangeItem();
                stormWaterInFlow.Quantity = qtdFlowFromStrom;
                stormWaterInFlow.ElementSet = stormWaterInflowNodes;
                inputExchangeItems.Add(stormWaterInFlow);
            }

            //
            //Properties
            //

            //Properties input / output exchange items
            for (int i = 1; i <= mohidLandEngine.GetNumberOfProperties(drainageNetworkInstanceID); i++)
            {
                int propertyIDNumber = mohidLandEngine.GetPropertyIDNumber(drainageNetworkInstanceID, i);
                string propertyName = mohidLandEngine.GetPropertyNameByIDNumber(propertyIDNumber);

                Quantity concentrationQuantity = new Quantity(unitConcentration, "Concentration of " + propertyName,
                                                              propertyIDNumber.ToString(), ValueType.Scalar, dimConcentration);

                OutputExchangeItem outletConc = new OutputExchangeItem();
                outletConc.Quantity = concentrationQuantity;
                outletConc.ElementSet = outletNode;

                outputExchangeItems.Add(outletConc);


                InputExchangeItem boundaryConc = new InputExchangeItem();
                boundaryConc.Quantity = concentrationQuantity;
                boundaryConc.ElementSet = outletNode;

                inputExchangeItems.Add(boundaryConc);

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

            Console.WriteLine("MOHID Land SetValues: " + setValuesWatch.ElapsedMilliseconds.ToString());
            Console.WriteLine("MOHID Land GetValues: " + getValuesWatch.ElapsedMilliseconds.ToString());
            Console.WriteLine("MOHID Land Perform  : " + performStepWatch.ElapsedMilliseconds.ToString());

        }

        public string GetComponentDescription()
        {
            return "Mohid Land Hydrological Model.";
        }

        public string GetComponentID()
        {
            return "Mohid Land Model ID";
        }

        public global::OpenMI.Standard.ITime GetCurrentTime()
        {
            DateTime currentTime = mohidLandEngine.GetCurrentTime();
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(currentTime));
        }

        public global::OpenMI.Standard.ITimeStamp GetEarliestNeededTime()
        {
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(mohidLandEngine.GetCurrentTime().AddSeconds(mohidLandEngine.GetCurrentTimeStep())));
        }

        public global::OpenMI.Standard.ITime GetInputTime(string QuantityID, string ElementSetID)
        {
            return new TimeStamp(CalendarConverter.Gregorian2ModifiedJulian(mohidLandEngine.GetCurrentTime().AddSeconds(mohidLandEngine.GetCurrentTimeStep())));
        }

        public double GetMissingValueDefinition()
        {
            return -99;
        }

        /// <summary>
        /// Get values for:
        ///     flow at outlet
        ///     concentrations at outlet
        ///     ponded water column
        /// </summary>
        public IValueSet GetValues(string QuantityID, string ElementSetID)
        {

            getValuesWatch.Start();

            double[] returnValues;

            if (QuantityID == qtdOutflow.ID)
            {
                returnValues = new double[1];
                returnValues[0] = mohidLandEngine.GetOutletFlow(drainageNetworkInstanceID);
            }
            else if (QuantityID == qtdWaterColumn.ID)
            {

                returnValues = new double[numberOfWaterPoints];
                mohidLandEngine.GetPondedWaterColumn(runoffInstanceID, numberOfWaterPoints, ref returnValues);

            }
            else if (QuantityID == qtdFlowToStorm.ID)
            {
                int numberOfOutflowNodes = mohidLandEngine.GetNumberOfStormWaterOutFlowNodes(drainageNetworkInstanceID);
                returnValues = new double[numberOfOutflowNodes];
                mohidLandEngine.GetStormWaterOutflow(drainageNetworkInstanceID, numberOfOutflowNodes, ref returnValues);
            }
            else //Concentration of properties
            {
                returnValues = new double[1];
                
                int propertyID = Convert.ToInt32(QuantityID);

                returnValues[0] = mohidLandEngine.GetOutletFlowConcentration(drainageNetworkInstanceID, propertyID);
            }

            ScalarSet values = new ScalarSet(returnValues);

            getValuesWatch.Stop();

            return values;
        }


        /// <summary>
        /// Set values for:
        ///     water level at outlet
        ///     discharges
        ///     property concentrations
        /// </summary>
        public void SetValues(string QuantityID, string ElementSetID, global::OpenMI.Standard.IValueSet values)
        {

            setValuesWatch.Start();

            if (QuantityID == qtdOutletLevel.ID)
            {
                double waterLevel = ((ScalarSet)values).data[0];
                mohidLandEngine.SetDownstreamWaterLevel(drainageNetworkInstanceID, waterLevel);
            }
            else if (QuantityID == qtdDischarges.ID)
            {
                ScalarSet flow = (ScalarSet)values;
                double[] flowValues = flow.data;
                mohidLandEngine.SetStormWaterModelFlow(runoffInstanceID, flow.Count, ref flowValues);
            }
            else if (QuantityID == qtdFlowFromStrom.ID)
            {
                ScalarSet flow = (ScalarSet) values;
                double[] flowvalues = flow.data;
                mohidLandEngine.SetStormWaterInflow(drainageNetworkInstanceID, flow.Count, ref flowvalues);
            }
            else //Concentration of properties
            {

                double concentration = ((ScalarSet)values).data[0];
                int propertyID = Convert.ToInt32(QuantityID);
                mohidLandEngine.SetDownStreamConcentration(drainageNetworkInstanceID, propertyID, concentration);
            }

            setValuesWatch.Stop();


        }


        public bool PerformTimeStep()
        {
            performStepWatch.Start();
            mohidLandEngine.PerformTimeStep();
            performStepWatch.Stop();
            return true;
        }


        #endregion
    }
}
