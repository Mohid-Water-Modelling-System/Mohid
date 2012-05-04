using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MOHID.OpenMI.MohidWater.Wrapper
{
    public class MohidWaterLinkableComponent : Oatc.OpenMI.Sdk.Wrapper.LinkableEngine
    {
        public MohidWaterLinkableComponent()
        {
            _engineApiAccess = new MohidWaterEngineWrapper();
        }

        protected override void SetEngineApiAccess()
        {
            _engineApiAccess = new MohidWaterEngineWrapper();
        }
    }
}
