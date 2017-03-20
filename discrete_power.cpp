//
// Compute discretized power spectrum 
//
// Ouput:
//     Column 1: k bin center
//     Column 2: P(k)
//     Column 3: number of independent modes
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>

#include "power.h"

using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  //
  // Command-line options (Boost program_options)
  //
  options_description opt("gaussian_rsd [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "power spectrum filename")
    ("boxsize", value<double>()->default_value(1000.0),   "boxsize")
    ("nc", value<size_t>()->default_value(64), "number of grids per dimension")
    ("k-unit", value<int>()->default_value(1), "0: fundamental frequency, 1: h/Mpc")
    ("k-min", value<double>()->default_value(0.0))
    ("k-max", value<double>()->default_value(1.0))
    ("dk", value<double>()->default_value(0.01), "k bin width")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt;
    return 0;
  }

  const string filename= vm["filename"].as<string>();
  const double boxsize= vm["boxsize"].as<double>();
  const int nc= vm["nc"].as<size_t>();

  const double knq= (M_PI*nc)/boxsize;
  const double knq_inv= boxsize/(M_PI*nc);


  const int k_unit= vm["k-unit"].as<int>();

  const double fac= 2.0*M_PI/boxsize;

  double k_min= vm["k-min"].as<double>();
  double k_max= vm["k-max"].as<double>();
  double dk=    vm["dk"].as<double>();
  
  if(k_unit == 0) {
    k_min *= fac;
    k_max *= fac;
    dk *= fac;
  }

  cerr << k_min << " " << k_max << " " << dk << endl;
  
  const int nbin= (int) round((k_max - k_min)/dk);


  PowerSpectrum* ps = power_alloc(filename.c_str());
  
  //
  // Allocate memory
  //
  int* const nmodes= (int*) calloc(sizeof(int), nbin);
  double* const kmean= (double*) calloc(sizeof(double), nbin);
  double* const P= (double*) calloc(sizeof(double), nbin);
  
  for(int ix=0; ix<nc; ++ix) {
   double kx= ix <= nc/2 ? fac*ix : fac*(ix - nc);

   for(int iy=0; iy<nc; ++iy) {
    double ky= iy <= nc/2 ? fac*iy : fac*(iy - nc);

    int iz0 = !(kx > 0.0f || (kx == 0.0 && ky > 0.0));

    // Avoid double counting on kz=plain
    // k=(0,0,0) dropped because this is 0
    // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
    //      1 otherwize
    //
      
    for(int iz=iz0; iz<nc/2+1; ++iz) {
      double kz= fac*iz;
      
      double k= sqrt(kx*kx + ky*ky + kz*kz);

      int i= (int) floor((k - k_min)/dk);
      
      if(0 <= i && i < nbin) {
	nmodes[i]++;
	kmean[i] += k;
	P[i] += power(ps, k);
      }
    }  
   }
  }

  for(int i=0; i<nbin; ++i) {
    if(nmodes[i] > 0) {
      printf("%e %e %d\n",
	     k_min + (i + 0.5)*dk,
	     P[i]/nmodes[i],
	     nmodes[i]);

    }
  }
  
  return 0;
}

