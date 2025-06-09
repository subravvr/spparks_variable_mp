/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_AM_THERMALFIELD_H
#define SPK_AM_THERMALFIELD_H

#include <math.h>

namespace RASTER_THERMAL {

namespace thermalfield {

   class ThermalField : public FieldShape {

      private:
         // Parameters
         // Const is good!
         const double _field_x, _field_y, _field_z, _ref_x, _ref_y, _ref_z;

      public:
         ThermalField(double fx, double fy, double fz, double rx, double ry, double rz):
         	_field_x(fx), _field_y(fy), _field_z(fz), _ref_x(rx), _ref_y(ry), _ref_z(rz) {}

         virtual ~ThermalField(){}


         bool is_inside(const double *XYZ) const {

					 if ((XYZ[0] < _ref_x) || (XYZ[0] > _ref_x + _field_x) || (XYZ[1] < _ref_y) || (XYZ[1] > _ref_y + _field_y) || (XYZ[2] < _ref_z) || (XYZ[2] > _ref_z + _field_z)) {
						 return false;
					 } else {
							return true;
						}

        }

        //This is probably the better place to do all the checks that we're currently doing in app_update

   };

}

}

#endif
