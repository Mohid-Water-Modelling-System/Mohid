using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Oatc.OpenMI.Sdk.Spatial;
using Oatc.OpenMI.Sdk.Wrapper;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.DevelopmentSupport;
using OpenMI.Standard;
using System.Collections;
using ValueType = OpenMI.Standard.ValueType;

namespace MOHID.OpenMI.MohidWater.Wrapper
{
    public class MohidWaterEngineWrapper : IEngine
    {

        #region Fields

        private MohidWaterEngineDotNetAccess MohidWaterEngine;
        private IList<InputExchangeItem> inputExchangeItems;
        private IList<OutputExchangeItem> outputExchangeItems;

        //TODO: Check how we can the correct instance IDs from the MOHID models. Implement a selector in the main, which receives the model ID?
        private int dischargeInstanceID = 1;
        private int horizontalGridInstanceID = 1;
        private int horizontalMapInstanceID = 1;
        private int hydrodynamicInstanceID = 1;
        private int waterPropertiesInstanceID = 1;

        private Quantity qtdDischargeFlow;      //Flow at select discharge points   -> input
        private Quantity qtdWaterLevel;         //Waterlevel at grid points         -> output

        private IList<Quantity> qtdProperties = new List<Quantity>();

        private int numberOfWaterPoints;

        private double[] modelGridValues1D;


        #endregion

        #region IEngine Members

        public void Initialize(System.Collections.Hashtable properties)
        {
            //List of exchange Items
            inputExchangeItems = new List<InputExchangeItem>();
            outputExchangeItems = new List<OutputExchangeItem>();

            //Initializes Engine
            MohidWaterEngine = new MohidWaterEngineDotNetAccess();
            MohidWaterEngine.Initialize(properties["FilePath"].ToString());

            //
            //Dimensions
            //
            Dimension dimFlow = new Dimension();
            dimFlow.SetPower(DimensionBase.Length, 3);
            dimFlow.SetPower(DimensionBase.Time, -1);

            Dimension dimWaterlevel = new Dimension();
            dimWaterlevel.SetPower(DimensionBase.Length, 1);

            Dimension dimConcentration = new Dimension();
            dimConcentration.SetPower(DimensionBase.Mass, 1);
            dimConcentration.SetPower(DimensionBase.Length, -3);

            //
            //Units
            //
            Unit unitFlow = new Unit("m3/sec", 1, 0, "cubic meter per second");
            Unit unitWaterLevel = new Unit("m", 1, 0, "sea water level");
            Unit unitConcentration = new Unit("mg/l", 1, 0, "miligram per liter");

            //
            //Quantities
            //
            qtdDischargeFlow = new Quantity(unitFlow, "Input Discharge Flow", "Discharge Flow",
                                            global::OpenMI.Standard.ValueType.Scalar, dimFlow);
            qtdWaterLevel = new Quantity(unitWaterLevel, "Waterlevel of the water surface", "Waterlevel",
                                         global::OpenMI.Standard.ValueType.Scalar, dimWaterlevel);

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
            for (int i = 1; i <= MohidWaterEngine.GetIUB(horizontalGridInstanceID); i++)
            {
                for (int j = 1; j <= MohidWaterEngine.GetJUB(horizontalGridInstanceID); j++)
                {
                    if (MohidWaterEngine.IsWaterPoint(horizontalMapInstanceID, i, j))
                    {
                        String name = "i=" + i.ToString() + "/j=" + j.ToString();

                        double[] xCoords = new double[5];
                        double[] yCoords = new double[5];
                        MohidWaterEngine.GetGridCellCoordinates(horizontalGridInstanceID, i, j, ref xCoords, ref yCoords);

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

            //allocates waterlevels1D
            modelGridValues1D = new double[numberOfWaterPoints];

            //Discharge Points
            ElementSet dischargePoints = new ElementSet("Discharge Points", "Discharge Points", ElementType.XYPoint,
                                                        spatialReference);

            //Flow input exchange to discharges configured as OpenMI Discharges
            for (int i = 1; i <= MohidWaterEngine.GetNumberOfDischarges(dischargeInstanceID); i++)
            {
                if (MohidWaterEngine.GetDischargeType(dischargeInstanceID, i) == 4)
                {
                    Element dischargeElement = new Element(MohidWaterEngine.GetDischargeName(dischargeInstanceID, i));
                    dischargeElement.AddVertex(
                        new Vertex(MohidWaterEngine.GetDischargeXCoordinate(dischargeInstanceID, i),
                                   MohidWaterEngine.GetDischargeYCoordinate(dischargeInstanceID, i), 0));
                    dischargePoints.AddElement(dischargeElement);
                }
            }

            //
            //Output Exchange Items
            //

            //Water Level of the Hydrodynamic model
            OutputExchangeItem waterlevel = new OutputExchangeItem();
            waterlevel.Quantity = qtdWaterLevel;
            waterlevel.ElementSet = modelGrid;
            outputExchangeItems.Add(waterlevel);


            //Properties of the Water properties model
            for (int idx = 1; idx <= MohidWaterEngine.GetNumberOfProperties(waterPropertiesInstanceID); idx++)
            {

                int propertyID = MohidWaterEngine.GetPropertyIDNumber(waterPropertiesInstanceID, idx);
                string propertyName = MohidWaterEngine.GetPropertyNameByIDNumber(propertyID);

                Quantity concentrationQuantity = new Quantity(unitConcentration, "Concentration of " + propertyName, propertyID.ToString(),
                                                              ValueType.Scalar, dimConcentration);

                qtdProperties.Add(concentrationQuantity);

                OutputExchangeItem concExchangeItem = new OutputExchangeItem();
                concExchangeItem.Quantity = concentrationQuantity;
                concExchangeItem.ElementSet = modelGrid;

                outputExchangeItems.Add(concExchangeItem);
            }


            //Flow input exchange to discharges configured as OpenMI Discharges
            for (int i = 1; i <= MohidWaterEngine.GetNumberOfDischarges(dischargeInstanceID); i++)
            {
                if (MohidWaterEngine.GetDischargeType(dischargeInstanceID, i) == 4)
                {

                    String pointName = MohidWaterEngine.GetDischargeName(dischargeInstanceID, i);

                    ElementSet dischargePoint = new ElementSet("Discharge Point of MOHID Water", i.ToString(), ElementType.XYPoint, spatialReference);

                    InputExchangeItem inputDischarge = new InputExchangeItem();
                    inputDischarge.Quantity = qtdDischargeFlow;

                    Element element = new Element("Point: " + pointName);
                    element.AddVertex(new Vertex(MohidWaterEngine.GetDischargeXCoordinate(dischargeInstanceID, i), MohidWaterEngine.GetDischargeYCoordinate(dischargeInstanceID, i), 0));
                    dischargePoint.AddElement(element);

                    inputDischarge.ElementSet = dischargePoint;

                    inputExchangeItems.Add(inputDischarge);


                    for (int idx = 1; idx <= MohidWaterEngine.GetNumberOfDischargeProperties(dischargeInstanceID, i); idx++)
                    {
                        int propertyID = MohidWaterEngine.GetDischargePropertyID(dischargeInstanceID, i, idx);

                        string propertyName = MohidWaterEngine.GetPropertyNameByIDNumber(propertyID);

                        Quantity concentrationQuantity = new Quantity(unitConcentration, propertyName, propertyID.ToString(), global::OpenMI.Standard.ValueType.Scalar, dimConcentration);

                        InputExchangeItem inputExchangeItem = new InputExchangeItem();
                        inputExchangeItem.ElementSet = dischargePoint;
                        inputExchangeItem.Quantity = concentrationQuantity;

                        inputExchangeItems.Add(inputExchangeItem);

                    }

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


        /// <summary>
        /// Get values for:
        ///     water level
        ///     properties
        /// </summary>
        public IValueSet GetValues(string QuantityID, string ElementSetID)
        {

            double[] returnValues;

            if (QuantityID == qtdWaterLevel.ID)
            {
                MohidWaterEngine.GetWaterLevel1D(hydrodynamicInstanceID, numberOfWaterPoints,  ref modelGridValues1D);

                return new ScalarSet(modelGridValues1D);
            }
            else //Gets concentration from Waterproperties
            {

                int propertyID = Convert.ToInt32(QuantityID);

                MohidWaterEngine.GetConcentration1D(waterPropertiesInstanceID, propertyID, numberOfWaterPoints, ref modelGridValues1D);

                return new ScalarSet(modelGridValues1D);

            }

        }


        public void SetValues(string QuantityID, string ElementSetID, IValueSet values)
        {
            if (QuantityID == qtdDischargeFlow.ID)
            {
                double flow = ((ScalarSet)values).data[0];
                MohidWaterEngine.SetDischargeFlow(dischargeInstanceID, Convert.ToInt32(ElementSetID), flow);
            }
            else //Set discharge concentrations
            {
                
                double concentration = ((ScalarSet)values).data[0];

                int dischargeID = Convert.ToInt32(ElementSetID);
                int propertyID = Convert.ToInt32(QuantityID);
                
                MohidWaterEngine.SetDischargeConcentration(dischargeInstanceID, dischargeID, propertyID, concentration);
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
