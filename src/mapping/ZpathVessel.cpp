/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class ZpathVessel : public vesselbase::FunctionVessel {
private:
  double invlambda;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  ZpathVessel( const vesselbase::VesselOptions& da );
  std::string function_description();
  bool calculate( std::vector<double>& buffer );
  double finalTransform( const double& val, double& dv );
};

PLUMED_REGISTER_VESSEL(ZpathVessel,"ZPATH")

void ZpathVessel::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords(keys);
}

void ZpathVessel::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("ZPATH",false,"calculate the distance from the low dimensionality manifold",true);
  keys.addOutputComponent("zpath","ZPATH","the distance from the path");
}

ZpathVessel::ZpathVessel( const vesselbase::VesselOptions& da ):
FunctionVessel(da)
{
  Mapping* mymap=dynamic_cast<Mapping*>( getAction() );
  plumed_massert( mymap, "ZpathVessel should only be used with mappings");
  invlambda = 1.0 / mymap->getLambda(); usetol=true;
}

std::string ZpathVessel::function_description(){
  return "the distance from the low-dimensional manifold";
}

bool ZpathVessel::calculate( std::vector<double>& buffer ){
  return addToBuffers( 1.0, 0.0, buffer );
}

double ZpathVessel::finalTransform( const double& val, double& dv ){
   dv = -invlambda / val;
   return -invlambda*std::log( val );
} 

}
}
