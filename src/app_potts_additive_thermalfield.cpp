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
#include "app_potts_additive_thermalfield.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "am_thermalfield.h"
#include "am_thermalfield_raster.h"
#include "error.h"

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <valarray>
#include <limits>

using namespace SPPARKS_NS;
using RASTER::Point;

/* ---------------------------------------------------------------------- */

AppPottsAdditiveThermalField::AppPottsAdditiveThermalField(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg),  MobilityOut(0)
{

   // only error check for this class, not derived classes
   if (strcmp(arg[0],"additive_thermalfield") == 0 && narg != 7 )
    error->all(FLERR,"Illegal app_style command");

   nspins = atoi(arg[1]); //Number of spins
   field_x = atof(arg[2]); //Length of the thermal field in x
   field_y = atof(arg[3]); //Length of the thermal field in y
   field_z = atof(arg[4]); //Length of the thermal field in z
	 exp_factor = atof(arg[5]); //Arrhenius prefactor value
	 activation_energy = atof(arg[6]); //Material specific activation energy
   filename = atof(arg[7]); //Name of the text file with temperature info

   //Define the layer object, this might work better in init_app
   ndouble = 1;
   allow_app_update = 1;
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAdditiveThermalField::input_app(char *command, int narg, char **arg)
{
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else error->all(FLERR,"Unrecognized command");
}


/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAdditiveThermalField::grow_app()
{
  spin = iarray[0];
  MobilityOut = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAdditiveThermalField::init_app()
{
   // Run base class init_app
   init_app_am();

   // Compute distance function based upon initial pool position
   app_update(0.0);

}

/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */

 /* ----------------------------------------------------------------------
    Read file with temperature values
 ------------------------------------------------------------------------- */
 std::valarray<double> AppPottsAdditiveThermalField::get_temperature_values(string filename, int field_x, int field_y, int field_z) {
 	using namespace std;
	int num_points = field_x * field_y * field_z;
 	valarray<double> temp_field(field_x*field_y*field_z);
	double current_t;
	int i = 0;

 	ifstream inFile(filename,ios::in);
 	if (!inFile.is_open()) {
		std:cerr << "There was a problem opening the temperature data file.\n";
		exit(1);
	}

 	while (inFile >> current_t)
	{
 		temp_field[i] = current_t;

		i++;
 	}
 	inFile.close();
 	return temp_field;

 }

void AppPottsAdditiveThermalField::app_update(double dt)
{
   // Move pool
   bool moved=app_update_am(dt);
   if(!moved)
      return;


	 Path current;
	 Point current_layer_dir = current.get_start();
	 string debug_string1("start:");
	 std::cout << debug_string1 << current_layer_dir << std::endl;
	 double origin[] = {0,0,0};
	 Point refpoint = compute_position_relative_to_pool(origin);
	 Point refpoint_rounded(std::round(refpoint[0]),std::round(refpoint[1]),std::round(refpoint[2]));


	 RASTER_THERMAL::thermalfield::ThermalField tf(field_x, field_y, field_z, -refpoint_rounded[0], -refpoint_rounded[1],-refpoint_rounded[2]);

	//Use the new position as input to the mobility calculation
	//Loop through all of the local sites and assign the new mobilities
	int num_points = field_x * field_y * field_z;
	// Get temperature from file
	std::ifstream inFile;
	std::valarray<double> temp_field = get_temperature_values("rosenthal_field_test.txt", field_x, field_y, field_z);


	for(int i=0;i<nlocal;i++){

		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};
		// Find index of point inside the thermal domain

		bool insidefield = tf.is_inside(XYZ);
		if (!insidefield){
			MobilityOut[i] = 0.0;
		} else{
			//string debug_string3("This site is inside the thermal domain: ");
			//std::cout << debug_string3 << i << std::endl;
			int thermal_index = 0;
			for(int thermal_z_increment=0; thermal_z_increment<field_z; thermal_z_increment++){
				for(int thermal_y_increment=0; thermal_y_increment<field_y; thermal_y_increment++){
					for(int thermal_x_increment=0; thermal_x_increment<field_x; thermal_x_increment++){
						if ( (XYZ[0]==thermal_x_increment-refpoint_rounded[0]) && (XYZ[1]==thermal_y_increment-refpoint_rounded[1]) && (XYZ[2]==thermal_z_increment-refpoint_rounded[2])){
							// Compute mobility
							MobilityOut[i] = compute_mobility(i,thermal_index,temp_field);
							//string debug_string2("Mobility at this site is ");
							//std::cout << debug_string2 << MobilityOut[i] << std::endl;
						}
						thermal_index+=1;
					}
				}
			}
		}





		}

	}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double AppPottsAdditiveThermalField::compute_mobility(int site, int thermal_index, std::valarray<double> temp_field)  {

	//Get the temperature at the current site in thermal domain
	double temp = temp_field[thermal_index];
	//Get the max temp in the domain
	double max_temp = temp_field.max();
	double site_mobility = exp_factor*exp(-activation_energy / (8.314 * temp));
	double max_mobility = exp_factor*exp(-activation_energy / (8.314 * max_temp));
	double normalized_mobility = site_mobility/max_mobility;
	//gradient in app_update, so here we'll just calculate the mobility based on distance

	MobilityOut[site] = site_mobility;

	return MobilityOut[site];
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
void AppPottsAdditiveThermalField::site_event_rejection(int site, RandomPark *random) {
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
				//string debug_string2("Reassigned spin to 1. Trying again: ");
        spin[site] = (int) (nspins*random->uniform());
				//std::cout << debug_string2 << spin[site] << std::endl;
      }
      return;
   }

   // events = spin flips to neighboring site different than self
   int j,m,value;
   int nevent = 0;
   int z = xyz[site][2];

   if((MobilityOut[site] > 0.0) && (MobilityOut[site] < 1.0) && (spin[site]==1)) {
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
