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

/* ----------------------------------------------------------------------
   Contributing author: Theron Rodgers and John Mitchell (Sandia)
------------------------------------------------------------------------- */

#include "string.h"
#include "math.h"
#include "app_potts_additive_spline.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "am_spline_mp.h"
#include "am_spline_raster.h"
#include "potts_am_path_parser.h"
#include "error.h"

#include <iostream>
#include <limits>
#include <valarray>
#include <fstream>
#include <string>

using namespace SPPARKS_NS;
using RASTER::Point;

/* ---------------------------------------------------------------------- */

AppPottsAdditiveSpline::AppPottsAdditiveSpline(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg),  MobilityOut(0)
{

   // only error check for this class, not derived classes
   if (strcmp(arg[0],"additive") == 0 && narg != 13 )
    error->all(FLERR,"Illegal app_style command");

   nspins = atoi(arg[1]); //Number of spins
   spot_width = atof(arg[2]); //Width of the melt pool
   melt_tail_length = atof(arg[3]); //Length of tail from meltpool midpoint
   melt_depth = atof(arg[4]); //How many lattice sites deep the melt pool is
   melt_cap = atof(arg[5]); //How many sites high the melt pool cap is
   cap_height = atof(arg[6]); //Height of the cap leading the meltpool
   HAZ = atof(arg[7]); //Size of the HAZ surrounding the melt pool (must be larger than spot_width)
   tail_HAZ = atof(arg[8]); //Length of hot zone behind meltpool (must be larger than melt_tail_length)
   depth_HAZ = atof(arg[9]); //Depth of the hot zone underneath the meltpool (must be larger than melt_depth)
   cap_HAZ = atof(arg[10]); //Size of HAZ infront of the melt pool (must be larger than cap_height)
   spline_coeffs_filename = atof(arg[11]); // filename containing spline coefficients
   knot_k1 = atof(arg[12]); // k1 parameter value related to the middle knot in the spline
   exp_factor = atof(arg[13]); //Exponential parameter for mobility decay in haz M(d) = exp(-exp_factor * d)

   //Define the layer object, this might work better in init_app
   ndouble = 1;
   allow_app_update = 1;
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAdditiveSpline::input_app(char *command, int narg, char **arg)
{
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAdditiveSpline::grow_app()
{
  spin = iarray[0];
  MobilityOut = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAdditiveSpline::init_app()
{
   // Run base class init_app
   init_app_am();
   // Compute distance function based upon initial pool position
   app_update(0.0);
}

std::valarray<double>  AppPottsAdditiveSpline::read_spline_vals(string spline_coeffs_filename) {
	using namespace std;
	valarray<double> spline_coeffs(8);
	double current_coeff;
	int i = 0;
	ifstream inFile1(spline_coeffs_filename,ios::in);
 	if (!inFile1.is_open()) {
		std:cerr << "There was a problem opening the depth data file.\n";
		exit(1);
	}
	while (inFile1 >> current_coeff)
	{
 		spline_coeffs[i] = current_coeff;
		i++;
 	}
	inFile1.close();
	return spline_coeffs;

}

/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */

void AppPottsAdditiveSpline::app_update(double dt)
{
   // read in spline coeffs
   std::ifstream inFile1;
   string spline_coeffs_filename = "concave.txt";
	std::valarray<double> spline_coeffs = read_spline_vals(spline_coeffs_filename);
   // Move pool
   bool moved=app_update_am(dt);
   if(!moved)
      return;


	//Use the new position as input to the mobility calculation
	//Loop through all of the local sites and assign the new mobilities

   // read in spline_coeffs

	//Specify the shape of the melt pool and then calculate the distance at each local site.
	RASTER_SPLINE::spline_pool::SplinePool sp(spot_width, melt_depth, spline_coeffs, knot_k1, melt_cap, melt_tail_length, cap_height, HAZ, tail_HAZ);

	//Go through all the local sites and calculate the distance.
   double d;
	for(int i=0;i<nlocal;i++){

		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};
		// Lattice point location relative to 'pool' position
		Point xyz_r_p= compute_position_relative_to_pool(XYZ);

		//Temporary assignment of xo, xo is in the melt pool's reference frame!
		double xo[]={xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};


		if (sp.is_inside(xo)) {
         MobilityOut[i] = 1;
      }
      else{
         MobilityOut[i] =0;
      }

		//Only calculate mobilities for things inside the HAZ bounds and below the active layer
	// 	if (d >= 0) {
	// 		MobilityOut[i] =  compute_mobility(i, d);
	// 	}
	// 	//Inside the pool so set mobilty to 1 (which randomizes things)
	// 	else if (d>=-5) {

	// 		MobilityOut[i] = 1;
	// 	}
	// 	//If we're outside the region of interest, make Mobility zero
	// 	else {
	// 		MobilityOut[i] = 0;
	// 	}
	// }
   }
}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double AppPottsAdditiveSpline::compute_mobility(int site, double d)  {

	//We're going to take care of categorizing all the little details of the mobility
	//gradient in app_update, so here we'll just calculate the mobility based on distance
	MobilityOut[site] = exp(-exp_factor * d);

	return MobilityOut[site];
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
void AppPottsAdditiveSpline::site_event_rejection(int site, RandomPark *random) {
   int oldstate = spin[site];
   double einitial = site_energy(site);

   if (MobilityOut[site] < 0.0) {
      MobilityOut[site] = 0.0;
      return;
   }

   if(MobilityOut[site] >= 1.0){
      //Mobility = 0.0;
      spin[site] = (int) (nspins*random->uniform());
      // While loop ensures spin[site] cannot be assigned 1
      while (spin[site] == 1) {
         spin[site] = (int) (nspins*random->uniform()); 
      }
      return;
   }

   // Ensures unmelted sites in HAZ stay unmelted
   if((MobilityOut[site] > 0.0) && (MobilityOut[site] < 1.0) && spin[site] == 1) {
      MobilityOut[site] = 0.0;
      return;
   }

   // events = spin flips to neighboring site different than self
   int j,m,value;
   int nevent = 0;
   int z = xyz[site][2];

   if((MobilityOut[site] > 0.0) && (MobilityOut[site] < 1.0)) {
      //(spin[i] != nspins) another criteria to exclude gg interaction
      for (j = 0; j < numneigh[site]; j++) {
         value = spin[neighbor[site][j]];
         if (value == spin[site] || value == nspins || value == 1) continue;
            for (m = 0; m < nevent; m++)
               if (value == unique[m]) break;
            if (m < nevent) continue;
         unique[nevent++] = value;
      }

      if (nevent == 0) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;
         spin[site] = unique[iran];
      double efinal = site_energy(site);

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > MobilityOut[site]){
            spin[site] = oldstate;
         }
      } else if (temperature == 0.0) {
         spin[site] = oldstate;
      } else if (random->uniform() > MobilityOut[site] * exp((einitial-efinal)*t_inverse)) {
         spin[site] = oldstate;
      }

      if (spin[site] != oldstate) naccept++;

      // set mask if site could not have changed
      // if site changed, unset mask of sites with affected propensity
      // OK to change mask of ghost sites since never used

      if (Lmask) {
         if (einitial < 0.5*numneigh[site]) mask[site] = 1;
         if (spin[site] != oldstate)
         for (int j = 0; j < numneigh[site]; j++)
            mask[neighbor[site][j]] = 0;
      }
   }
}
