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

#ifndef SPK_AM_SPLINE_H
#define SPK_AM_SPLINE_H

#include <math.h>
#include <valarray>

namespace RASTER_SPLINE {

namespace spline_pool {

   class SplinePool : public SplinePoolShape {

      private:
         // Parameters 
         // Const is good!
         const double _p_width, _p_depth, _p_cap, _tail_length, _cap_height, _HAZ_width, _tail_HAZ,_knot_k1;
		 std::valarray<double> _spline_coeffs;
		 

      public:
         SplinePool(double p_width, double p_depth, std::valarray<double> spline_coeffs, double knot_k1, double p_cap, double tail_length, double cap_height, double HAZ_width, double tail_HAZ):
         	_p_width(p_width), _p_depth(p_depth), _spline_coeffs(spline_coeffs), _p_cap(p_cap), _tail_length(tail_length), _cap_height(cap_height), _HAZ_width(HAZ_width), _tail_HAZ(tail_HAZ), _knot_k1(knot_k1) {}

         virtual ~SplinePool(){}
         
         
         bool is_inside(const double *XYZ) const {
        	double active_p_length;
	  		double active_p_height;
        	
        	//If z is positive, we're above the current layer and can't be inside
			if (XYZ[2] > 0) {
			 	active_p_height = _p_cap;
				if (XYZ[0] < 0) {
			 	active_p_length = _tail_length;	
			 	}
			 	else {
			 	active_p_length = _cap_height;
			 	}
			 return pow(XYZ[1]/(_p_width * 0.5),2) + pow(XYZ[0]/active_p_length,2) + pow(XYZ[2]/active_p_height,2) < 1.0;
			 }

			 else{
				if (XYZ[0] < 0) {
			 	active_p_length = _tail_length;	
			 	}
			 	else {
			 	active_p_length = _cap_height;
				}

				// check spline depth here
				double a1 = _spline_coeffs[0];
				double b1 = _spline_coeffs[1];
				double c1 = _spline_coeffs[2];
				double d1 = _spline_coeffs[3];
				double a2 = _spline_coeffs[4];
				double b2 = _spline_coeffs[5];
				double c2 = _spline_coeffs[6];
				double d2 = _spline_coeffs[7];

				// double effective_depth = _p_depth*sqrt(1-pow(XYZ[0]/active_p_length,2));
				// std::cout << effective_depth << std::endl;
				// double scaling_factor = effective_depth/_p_depth;
				double scaling_factor = 1;

				double yval = abs(XYZ[1]); // to account for symmetry
				if (yval<_knot_k1*_p_width/2){
					// evaluate spline 1
					double spline_depth = a1*pow(yval,3) + b1*pow(yval,2) + c1*yval + d1;
					if (XYZ[2]>= spline_depth*scaling_factor){
						return true;
					}
					else{
						return false;
					}
				}
				else{
					// evaluate spline 2
					double spline_depth = a2*pow(yval,3) + b2*pow(yval,2) + c2*yval + d2;
					if (XYZ[2]>= spline_depth*scaling_factor){
						return true;
					}
					else{
						return false;
					}
				}
			 }
            
            
         }


         virtual double distance(const double *XYZ) const {

			 double active_p_length;
		   	 double active_p_height;
			 int vector_length;
			 double d = -1.0;
			 double working_tail_haz = _tail_HAZ - _tail_length;

			 			 
			 //Test to see if we're inside the melt pool, return -1 if so
			 if(is_inside(XYZ)) return d;
			 
			 //Change active_p_length depending on whether XYZ[1] is pos or neg
			 if (XYZ[0] < 0) {
			 	active_p_length = _tail_length;	
			 }
			 else {
			 	active_p_length = _cap_height;
			 }
			 
			 //Change active_p_height depending on whether XYZ[2] is pos or neg
			 if (XYZ[2] < 0) {
			 	active_p_height = _p_depth;	
			 }
			 else {
			 	active_p_height = _p_cap;
			 }			
			 
			 //Determine the length from the melt pool center to our point.
			 //This will be our maximum possible iteration
			 vector_length = sqrt(XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1] + XYZ[2] * XYZ[2]);
			 	

			 
			//Calculate distance from trailing ellipsoid.
			//We should be able to return distance for any coordinate (although this would be expensive for far away points)
			for(int i = 0; i < vector_length; i++) {
			
				if (pow(XYZ[1]/(_p_width * 0.5 + i),2) + pow(XYZ[0]/(active_p_length + i),2) + pow(XYZ[2]/(active_p_height + i),2) <= 1) {
					
					//Smallest shell that the point is inside
					d = i;
					break;
				}
			}
							 	
            return d;
        }
        
        //This is probably the better place to do all the checks that we're currently doing in app_update
        
   };

}

}

#endif
