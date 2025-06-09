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

#ifdef APP_CLASS
AppStyle(potts/additive_spline,AppPottsAdditiveSpline)

#else

#ifndef SPK_APP_POTTS_ADDITIVE_SPLINE
#define SPK_APP_POTTS_ADDITIVE_SPLINE

#include <stdlib.h>
#include <vector>
#include <valarray>
#include <string>
#include "app_potts.h"
#include "am_spline_raster.h"
#include "potts_am_path_parser.h"

//using std::map;
//using RASTER::Pass;
//using RASTER::Path;
//using RASTER::Layer;

namespace SPPARKS_NS {

class AppPottsAdditiveSpline : public PottsAmPathParser {
 public:
  AppPottsAdditiveSpline(class SPPARKS *, int, char **);
  virtual void grow_app();
  virtual void init_app();
  virtual void site_event_rejection(int, RandomPark *);
  void input_app(char *, int , char **);
  double compute_mobility(int, double);
  void app_update(double);
  std::valarray<double> read_spline_vals(string);

 //Remove all of the variables we don't actually need
 protected:
 double *MobilityOut;

 double spot_width;
 string spline_coeffs_filename;
 double knot_k1;
 double melt_tail_length;
 double tail_HAZ;
 double exp_factor;
 double cap_height;
 double HAZ;
 double melt_depth;
 double melt_cap;
 double cap_HAZ;
 double depth_HAZ;
// private:
//   double build_layer_z;
//   map<int,Pass> passes;
//   map<int,Path> paths;
//   std::vector<Layer> layers;
//   std::vector<Layer>::iterator active_layer;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
