/**************************************************************************** 

  Copyright (C) 2011 Atsuto Seko
  seko@cms.mtl.kyoto-u.ac.jp

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to
  the Free Software Foundation, Inc., 51 Franklin Street,
  Fifth Floor, Boston, MA 02110-1301, USA, or see
    http://www.gnu.org/copyleft/gpl.txt

    Main program for calculating electrostatic energy using Ewald method

******************************************************************************/

#include <iostream>
#include <vector>

#include "input.h"
#include "output.h"
#include "ewald.h"
#include "include/math.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

int main(int argc, char *argv[]){

    // reading input files
    Input ip;
    ip.input_ewald(argc, argv);

    dmatrix axis = ip.get_axis();
    std::vector<dvector> position = ip.get_position();

    std::vector<double> charge = ip.get_charge();

    double volume = ip.get_volume();
    dmatrix reciprocal_axis = ip.get_reciprocal_axis();

    double eta = ip.get_eta();
    double rmax = ip.get_rmax();
    double gmax = ip.get_gmax();

    Ewald ewald
        (axis, position, volume, reciprocal_axis, charge, eta, rmax, gmax);
    ewald.calc_energy();

    double pot_r = ewald.get_real_energy();
    double pot_g = ewald.get_reciprocal_energy();
    double pot_s = ewald.get_self_energy();
    double energy = ewald.get_energy();
    std::cout.precision(18);
    std::cout << "  energy (real space) = " << pot_r << std::endl;
    std::cout << "  energy (reciprocal space) = " << pot_g << std::endl;
    std::cout << "  energy (self) = " << pot_s << std::endl;
    std::cout << "  electrostatic energy = " << energy << std::endl;

    return 0;

}
