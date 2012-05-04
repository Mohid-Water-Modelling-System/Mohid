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
using OpenMI.Standard;
using Oatc.OpenMI.Sdk.Backbone;

namespace MOHID.OpenMI.UnitTest
{
    /// <summary>
    /// Summary description for LinkFactory.
    /// </summary>
    public class LinkFactory
    {
        public LinkFactory()
        {
            //
            // TODO: Add constructor logic here
            //
        }

        public static ILink CreateLink(ILinkableComponent sourceComponent, string sourceQuantityID, string sourceElementSetID, ILinkableComponent targetComponent, string targetQuantityID, string targetElementSetID, string[] dataOperationIDs)
        {
            string linkID = sourceComponent.ComponentID + "(" + sourceQuantityID + ", " + sourceElementSetID + ") to " + targetComponent.ComponentID + "(" + targetQuantityID + ", " + targetElementSetID + ")";

            int outputExchangeItemIndex = -1;
            for (int i = 0; i < sourceComponent.OutputExchangeItemCount; i++)
            {
                if (sourceComponent.GetOutputExchangeItem(i).Quantity.ID == sourceQuantityID && sourceComponent.GetOutputExchangeItem(i).ElementSet.ID == sourceElementSetID)
                {
                    outputExchangeItemIndex = i;
                }
            }

            if (outputExchangeItemIndex < 0)
            {
                throw new Exception("Exception in LinkFactory.CreateLink, failed to find output exchange item for " + sourceQuantityID + ", " + sourceElementSetID + " during creation of link: " + linkID);
            }

            int inputExchangeItemIndex = -1;
            for (int i = 0; i < targetComponent.InputExchangeItemCount; i++)
            {
                if (targetComponent.GetInputExchangeItem(i).Quantity.ID == targetQuantityID && targetComponent.GetInputExchangeItem(i).ElementSet.ID == targetElementSetID)
                {
                    inputExchangeItemIndex = i;
                }
            }

            if (inputExchangeItemIndex < 0)
            {
                throw new Exception("Exception in LinkFactory.CreateLink, failed to find input exchange item for " + targetQuantityID + ", " + targetElementSetID + " during creation of link: " + linkID);
            }


            Link link = new Link();
            link.ID = linkID;
            link.Description = linkID;
            link.SourceComponent = sourceComponent;
            link.SourceQuantity = sourceComponent.GetOutputExchangeItem(outputExchangeItemIndex).Quantity;
            link.SourceElementSet = sourceComponent.GetOutputExchangeItem(outputExchangeItemIndex).ElementSet;
            link.TargetComponent = targetComponent;
            link.TargetQuantity = targetComponent.GetInputExchangeItem(inputExchangeItemIndex).Quantity;
            link.TargetElementSet = targetComponent.GetInputExchangeItem(inputExchangeItemIndex).ElementSet;

            bool dataOperationWasFound = false;

            for (int n = 0; n < dataOperationIDs.Length; n++)
            {
                dataOperationWasFound = false;
                for (int i = 0; i < sourceComponent.GetOutputExchangeItem(outputExchangeItemIndex).DataOperationCount; i++)
                {
                    if (sourceComponent.GetOutputExchangeItem(outputExchangeItemIndex).GetDataOperation(i).ID == dataOperationIDs[n])
                    {
                        link.AddDataOperation(sourceComponent.GetOutputExchangeItem(outputExchangeItemIndex).GetDataOperation(i));
                        dataOperationWasFound = true;
                    }
                }
                if (!dataOperationWasFound)
                {
                    throw new Exception("failed to find dataOperation: " + dataOperationIDs[n] + " during creation of link: " + linkID);
                }
            }

            return (ILink)link;
        }
    }
}
