using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MOHID.OpenMI.MohidLand.Wrapper
{
    public class MohidLandLinkableComponent : Oatc.OpenMI.Sdk.Wrapper.LinkableEngine
    {
        public MohidLandLinkableComponent()
        {
            _engineApiAccess = new MohidLandEngineWrapper();
        }

        protected override void SetEngineApiAccess()
        {
            _engineApiAccess = new MohidLandEngineWrapper();
        }


    }
}
