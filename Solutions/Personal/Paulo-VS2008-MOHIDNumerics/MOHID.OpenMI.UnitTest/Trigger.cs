#region Copyright
///////////////////////////////////////////////////////////
//
//    Copyright (C) 2005 HarmonIT Consortium
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    or look at URL www.gnu.org/licenses/lgpl.html
//
//    Contact info: 
//      URL: www.openmi.org
//	Email: sourcecode@openmi.org
//	Discussion forum available at www.sourceforge.net
//
//      Coordinator: Roger Moore, CEH Wallingford, Wallingford, Oxon, UK
//
///////////////////////////////////////////////////////////
//
//  Original author: Jan B. Gregersen, DHI - Water & Environment, Horsholm, Denmark
//  Created on:      December 1st 2005
//  Version:         1.0.0 
//
//  Modification history:
//  
//
///////////////////////////////////////////////////////////
#endregion

using System;
using System.Collections;
using OpenMI.Standard;
using Oatc.OpenMI.Sdk.Backbone;
using Oatc.OpenMI.Sdk.Buffer;

namespace MOHID.OpenMI.UnitTest
{
    /// <summary>
    /// Summary description for Trigger.
    /// </summary>
    public class Trigger : ILinkableComponent
    {
        private ILink _link;
        private Oatc.OpenMI.Sdk.Buffer.SmartBuffer _resultsBuffer;
        private TimeStamp _earliestInputTime;

        public Trigger()
        {
            _resultsBuffer = new SmartBuffer();
            _earliestInputTime = new TimeStamp(0);
        }

        public SmartBuffer ResultsBuffer
        {
            get
            {
                return _resultsBuffer;
            }
        }

        public void Finish()
        {
        }

        public ITimeStamp EarliestInputTime
        {
            get
            {
                return _earliestInputTime;
            }
        }

        public void AddLink(ILink link)
        {
            _link = link;
        }

        public void Dispose()
        {
            // TODO:  Add Trigger.Dispose implementation
        }

        public IValueSet GetValues(ITime time, string linkID)
        {
            // TODO:  Add Trigger.GetValues implementation
            return null;
        }

        public string ComponentDescription
        {
            get
            {
                // TODO:  Add Trigger.Description getter implementation
                return null;
            }
        }


        public string ComponentID
        {
            get
            {
                // TODO:  Add Trigger.ID getter implementation
                return null;
            }
        }


        public void Initialize(IArgument[] properties)
        {
            // TODO:  Add Trigger.Initialize implementation
        }

        public string ModelID
        {
            get
            {
                // TODO:  Add Trigger.ID getter implementation
                return "Trigger";
            }
        }

        public string ModelDescription
        {
            get
            {
                // TODO:  Add Trigger.ID getter implementation
                return null;
            }
        }


        public void Prepare()
        {
            // TODO:  Add Trigger.PrepareForComputation implementation
        }

        public void RemoveLink(string linkID)
        {
            // TODO:  Add Trigger.RemoveLink implementation
        }

        public int InputExchangeItemCount
        {
            get
            {
                return 1;
            }
        }

        public int OutputExchangeItemCount
        {
            get
            {
                return 0;
            }
        }

        public IInputExchangeItem GetInputExchangeItem(int index)
        {

            // -- create a flow quanitity --
            Dimension flowDimension = new Dimension();
            flowDimension.SetPower(DimensionBase.Length, 3);
            flowDimension.SetPower(DimensionBase.Time, -1);
            Unit literPrSecUnit = new Unit("LiterPrSecond", 0.001, 0, "Liters pr Second");
            Quantity flowQuantity = new Quantity(literPrSecUnit, "Flow", "Flow", global::OpenMI.Standard.ValueType.Scalar, flowDimension);

            Element element = new Element();
            element.ID = "DummyElement";
            ElementSet elementSet = new ElementSet("Dummy ElementSet", "DummyElementSet", ElementType.IDBased, new SpatialReference("no reference"));
            elementSet.AddElement(element);

            InputExchangeItem inputExchangeItem = new InputExchangeItem();
            inputExchangeItem.ElementSet = elementSet;
            inputExchangeItem.Quantity = flowQuantity;

            return inputExchangeItem;
        }


        public IOutputExchangeItem GetOutputExchangeItem(int index)
        {

            return null;


        }

        public ITimeSpan TimeHorizon
        {
            get
            {
                return null;
            }
        }

        public string Validate()
        {
            //TODO: Inplement this method correctly
            return "";
        }

        public void Run(ITime[] GetValuesTimes)
        {
            for (int i = 0; i < GetValuesTimes.Length; i++)
            {
                _resultsBuffer.AddValues(GetValuesTimes[i], (IScalarSet)_link.SourceComponent.GetValues(GetValuesTimes[i], _link.ID));
                _earliestInputTime.ModifiedJulianDay = ((ITimeStamp)GetValuesTimes[i]).ModifiedJulianDay;
            }
        }

        public void Run(TimeStamp runToTime)
        {
            //IScalarSet scalarSet = new ScalarSet();

            _earliestInputTime = runToTime;
            _link.SourceComponent.GetValues(runToTime, _link.ID);


            //ScalarSet scalarSet = new ScalarSet((IScalarSet)_link.SourceComponent.GetValues(time, _link.ID));
            //_earliestInputTime.ModifiedJulianDay = time.ModifiedJulianDay;
            //_resultsBuffer.AddValues(time, scalarSet);
        }

        public int GetResultsCount()
        {
            return _resultsBuffer.ValuesCount;
        }

        public double GetResult(int index)
        {
            return ((ScalarSet)_resultsBuffer.GetValuesAt(_resultsBuffer.TimesCount - 1)).GetScalar(index);
        }


        #region IPublisher Members

        public void SendEvent(IEvent Event)
        {
            // TODO:  Add Trigger.SendEvent implementation
        }

        public void UnSubscribe(IListener listener, EventType eventType)
        {
            // TODO:  Add Trigger.UnSubscribe implementation
        }

        public EventType GetPublishedEventType(int providedEventTypeIndex)
        {
            // TODO:  Add Trigger.GetPublishedEventType implementation
            return EventType.Informative;
        }

        public void Subscribe(IListener listener, EventType eventType)
        {
            // TODO:  Add Trigger.Subscribe implementation
        }

        public int GetPublishedEventTypeCount()
        {
            // TODO:  Add Trigger.GetPublishedEventTypeCount implementation
            return 0;
        }

        #endregion
    }
}
