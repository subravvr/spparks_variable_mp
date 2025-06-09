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
#include "am_artsurf_raster.h"
#include "math.h"
#include "app_potts_additive_artsurf.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "am_ellipsoid.h"
#include "potts_am_path_parser.h"
#include "error.h"

#include <iostream>
#include <limits>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <valarray>
#include <sstream>
#include <cmath>

using namespace SPPARKS_NS;
using RASTER::Point;

/* ---------------------------------------------------------------------- */

AppPottsAdditiveArtSurf::AppPottsAdditiveArtSurf(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg),  MobilityOut(0)
{

   // only error check for this class, not derived classes
   if (strcmp(arg[0],"additive") == 0 && narg != 17 )
    error->all(FLERR,"Illegal app_style command");

   nspins = atoi(arg[1]); //Number of spins
   melt_tail_length = atof(arg[2]); //Length of tail from meltpool midpoint
   melt_depth_filename = arg[3]; // depth sequence
   melt_width_filename = arg[4]; // width sequence
   melt_cap_filename = arg[5]; // cap height from scime
   hatch_val = atof(arg[6]);
   thickness_val = atof(arg[7]);
   cap_height = atof(arg[8]); 
   width_HAZ_multiplier = atof(arg[9]); // HAZ-W/W
   depth_HAZ_multiplier = atof(arg[10]); // HAD-D/D
   tail_HAZ = atof(arg[11]); //Length of hot zone behind meltpool (must be larger than melt_tail_length)
   cap_HAZ = atof(arg[12]); //Size of HAZ infront of the melt pool (must be larger than cap_height)
   exp_factor = atof(arg[13]); //Exponential parameter for mobility decay in haz M(d) = exp(-exp_factor * d)
   depth_length = atof(arg[14]); // length of depth sequence for memory allocation
   width_length = atof(arg[15]); // length of width sequence for memory allocation
   cap_length = atof(arg[16]); // length of cap sequence for memory allocation

   ticker = atof(arg[17]); // ticker for app update
   ticker_mod = atof(arg[18]); // mod for ticker

   //Define the layer object, this might work better in init_app
   ndouble = 1;
   allow_app_update = 1;
   recreate_arrays();

   depth_sequence = read_depth_values(melt_depth_filename);
   width_sequence = read_width_values(melt_width_filename);
   cap_sequence = read_cap_values(melt_cap_filename);
   
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAdditiveArtSurf::input_app(char *command, int narg, char **arg)
{
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAdditiveArtSurf::grow_app()
{
  spin = iarray[0];
  MobilityOut = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAdditiveArtSurf::init_app()
{
   // Run base class init_app
   init_app_am();

   // Compute distance function based upon initial pool position
   app_update(0.0);
}


std::valarray<double>  AppPottsAdditiveArtSurf::read_depth_values(string melt_depth_filename) {
	using namespace std;
	valarray<double> depth_sequence(depth_length);
	double current_d;
	int i = 0;
	ifstream inFile1(melt_depth_filename,ios::in);
 	if (!inFile1.is_open()) {
		std:cerr << "There was a problem opening the depth data file.\n";
		exit(1);
	}
	while (inFile1 >> current_d)
	{
 		depth_sequence[i] = std::round(current_d);
		i++;
 	}
	inFile1.close();
	return depth_sequence;

}
std::valarray<double>  AppPottsAdditiveArtSurf::read_width_values(string melt_width_filename) {
	using namespace std;
	valarray<double> width_sequence(width_length);
	double current_w;
	int j = 0;
	ifstream inFile2(melt_width_filename,ios::in);
 	if (!inFile2.is_open()) {
		std:cerr << "There was a problem opening the width data file.\n";
		exit(1);
	}
	while (inFile2 >> current_w)
	{
 		width_sequence[j] = std::round(current_w);
		j++;
 	}
	inFile2.close();
	return width_sequence;

}
std::valarray<double>  AppPottsAdditiveArtSurf::read_cap_values(string melt_cap_filename) {
	using namespace std;
	valarray<double> cap_sequence(cap_length);
	double current_c;
	int k = 0;
	ifstream inFile3(melt_cap_filename,ios::in);
 	if (!inFile3.is_open()) {
		std:cerr << "There was a problem opening the cap data file.\n";
		exit(1);
	}
	while (inFile3 >> current_c)
	{
 		cap_sequence[k] = std::round(current_c);
		k++;
 	}
	inFile3.close();
	return cap_sequence;

}




/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */

void AppPottsAdditiveArtSurf::app_update(double dt)
{  
	// initialize melt pool variables
	double melt_depth;
	double spot_width;
	double HAZ;
	double depth_HAZ;
	double melt_cap;

	int mod_tick = std::round(ticker/(ticker_mod));
	bool increment_flag = false;
	

	if (mod_tick==0){
		melt_depth = depth_sequence[0];
		spot_width = width_sequence[0];
		melt_cap = cap_sequence[0];
		// std::cout << melt_depth << std::endl;
		HAZ = width_HAZ_multiplier*spot_width;
		depth_HAZ = depth_HAZ_multiplier*melt_depth;
		//melt_cap = (4*hatch_val*thickness_val*packing_fraction - M_PI*spot_width*melt_depth)/(M_PI*spot_width);
	}


   
	//Use the new position as input to the mobility calculation
	//Loop through all of the local sites and assign the new mobilities

	// std::cout<<ticker<<std::endl;


	while ( std::abs(mod_tick - dt*ticker) == 3/26 && increment_flag==false){
		if (mod_tick>depth_length-1)
		{
		mod_tick = depth_length-1;
		melt_depth = depth_sequence[mod_tick];
    	spot_width = width_sequence[mod_tick];
		melt_cap = cap_sequence[mod_tick];
		// std::cout << "end of mp geometry info" << std::endl;
		HAZ = width_HAZ_multiplier*spot_width;
		depth_HAZ = depth_HAZ_multiplier*melt_depth;
		//melt_cap = (4*hatch_val*thickness_val*packing_fraction - M_PI*spot_width*melt_depth)/(M_PI*spot_width);
		}
		else{
		increment_flag = true;
		//std::cout << "perform increment" << std::endl;
		melt_depth = depth_sequence[mod_tick];
    	spot_width = width_sequence[mod_tick];
		melt_cap = cap_sequence[mod_tick];
		// std::cout << mod_tick << std::endl;
		HAZ = width_HAZ_multiplier*spot_width;
		depth_HAZ = depth_HAZ_multiplier*melt_depth;
		// melt_cap = (4*hatch_val*thickness_val*packing_fraction - M_PI*spot_width*melt_depth)/(M_PI*spot_width);
		}

		

		double origin[] = {0,0,0};
		Point pos = compute_position_relative_to_pool(origin);
		double x_coord_pool = -pos[0];
		//string debug_string1("pos : ");
		//string debug_string2(" depth : ");
		//std::cout << debug_string1 << x_coord_pool << debug_string2 << melt_depth << std::endl;
		//std::cout << melt_depth << std::endl;
	}
	if (std::abs(mod_tick-dt*ticker) != 2/26 ){
		increment_flag = false;
	}
	
	
	// Move pool
   bool moved=app_update_am(dt);
   if(!moved){
	return;
   }
   else{
	ticker+=1;
	mod_tick = std::round(ticker/ticker_mod);
   }




	//Specify the shape of the melt pool and then calculate the distance at each local site.
	RASTER::pool_shape::AmEllipsoid ae(spot_width, melt_depth, melt_cap, melt_tail_length, cap_height, HAZ, tail_HAZ);

	//Go through all the local sites and calculate the distance.
   double d;
	for(int i=0;i<nlocal;i++){

		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};
		// Lattice point location relative to 'pool' position
		Point xyz_r_p=compute_position_relative_to_pool(XYZ);

		//Temporary assignment of xo, xo is in the melt pool's reference frame!
		double xo[]={xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};


		if(xo[0] < 0 && ((xo[2] <= 0 && abs(xo[2]) <= melt_depth)||(xo[2]>0 && abs(xo[2]) <= melt_cap))  && xo[0] > -tail_HAZ && abs(xo[1]) <= spot_width/2.0) {

			//If we're in the fusion zone, calculate distance
			if (abs(xo[1]) <= spot_width * 0.5 && abs(xo[0]) <= tail_HAZ) {
				d = ae.distance(xo);
			}
			//If we're in the HAZ, calculate distance
			//else if (abs(xo[1]) <= HAZ/2.0 && abs(xo[0]) <tail_HAZ) {
				//d = ae.distance(xo);
			//}
		}
		//If we're in front of the pool, look out to a distance cap_HAZ away
		else if (abs(xo[0]) <= cap_HAZ && ((xo[2] <=0 && abs(xo[2]) <= melt_depth)||(xo[2]>0 && abs(xo[2]) <= melt_cap)) && abs(xo[1]) <= spot_width/2.0) {
			d = ae.distance(xo);
		}
		else d = -10;

		//Only calculate mobilities for things inside the HAZ bounds and below the active layer
		if (d >= 0) {
			MobilityOut[i] =  compute_mobility(i, d);
		}
		//Inside the pool so set mobilty to 1 (which randomizes things)
		else if (d > -5) {

			MobilityOut[i] = 1;
		}
		//If we're outside the region of interest, make Mobility zero
		else {
			MobilityOut[i] = 0;
		}
	}
}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double AppPottsAdditiveArtSurf::compute_mobility(int site, double d)  {

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
void AppPottsAdditiveArtSurf::site_event_rejection(int site, RandomPark *random) {
   int oldstate = spin[site];
   double einitial = 0;

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
   }

      if (nevent == 0) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;
         spin[site] = unique[iran];
      double efinal = 1;

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > MobilityOut[site]){
            spin[site] = oldstate;
         }
		 else if (temperature == 0.0) {
         spin[site] = oldstate;
		 } 
		 else if (random->uniform() > MobilityOut[site] * exp((einitial-efinal)*t_inverse)) {
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
