/**************************************************************************** 

  Copyright (C) 2013 Atsuto Seko
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

    Main program for Monte Carlo simulation

******************************************************************************/

#include <iostream>
#include <vector>
#include <algorithm>

#include "input.h"
#include "output.h"
#include "mc_initialize.h"
#include "mc_binary.h"
#include "gcmc_binary.h"
#include "mc.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;

int main(){

    Input ip;
    ip.input_mc();

    std::vector<int> n_type = ip.get_n_type();
    int n_sublattice = n_type.size();
    const std::string simulation_type = ip.get_simulation_type();
    const int i_ewald = ip.get_i_ewald();
    if (simulation_type == "gcmc" and n_sublattice > 1){
        std::cout << "gcmc is not available for multilattice CE." << std::endl;
        exit(8);
    }
    if (simulation_type == "gcmc" and i_ewald == 1) {
        std::cout << "Point-charge model is available for cmc." << std::endl;
        exit(8);
    }
    if (simulation_type != "cmc" and simulation_type != "gcmc"){
        std::cerr << " simulation_type must be cmc or gcmc." << std::endl;
        exit(8);
    }

    const MC_initialize mc_init(ip);

    // checking whether system is binary or not
    int binary(0);
    if (n_type.size() == count(n_type.begin(), n_type.end(), 2)) binary = 1;

    double energy, energy_ave, specific_heat;
    dvector correlation_ave;
    std::vector<int> spin;

    std::vector<std::vector<int> > n_change_array = ip.get_n_change_array();

    if (binary == 1 and simulation_type == "cmc"){
        MC_binary mc(ip, mc_init);
        if (i_ewald == 0) mc.cmc(mc_init);
        else if (i_ewald == 1) mc.cmc_ewald(ip, mc_init);
        energy = mc.get_energy();
        energy_ave = mc.get_energy_average();
        correlation_ave = mc.get_correlation_average();
        specific_heat = mc.get_specific_heat();
        spin = mc.get_spin();
    }
    else if (binary == 1 and simulation_type == "gcmc"){
        GCMC_binary mc(ip, mc_init);
        if (n_change_array.size() == 0) mc.gcmc(ip, mc_init);
        else if (n_change_array.size() > 1) mc.gcmc_pseudo(ip, mc_init);
        energy = mc.get_energy();
        energy_ave = mc.get_energy_average();
        correlation_ave = mc.get_correlation_average();
        specific_heat = mc.get_specific_heat();
        spin = mc.get_spin();
    }
    else if (binary == 0){
        MC mc(ip, mc_init);
        if (simulation_type == "cmc" and i_ewald == 0) mc.cmc(mc_init);
        else if (simulation_type == "gcmc" and i_ewald == 0) 
            mc.gcmc(ip, mc_init);
        else if (simulation_type == "cmc" and i_ewald == 1) {
            mc.cmc_ewald(ip, mc_init);
        }
        energy = mc.get_energy();
        energy_ave = mc.get_energy_average();
        correlation_ave = mc.get_correlation_average();
        specific_heat = mc.get_specific_heat();
        spin = mc.get_spin();
    }

    const std::vector<double> temp_array = ip.get_temp_array();
    const dmatrix axis = mc_init.get_axis();
    const std::vector<int> n_atoms = mc_init.get_n_atoms();
    const std::vector<dvector> position = mc_init.get_position();
    const std::vector<int> cluster_index = mc_init.get_cluster_index();
    std::vector<int> spin_array = ip.get_spin_array();

    Output op;
    op.output_mc
    (spin, energy, energy_ave, specific_heat, correlation_ave, 
     temp_array, axis, n_atoms, position, cluster_index, 
     n_type, spin_array);

    return 0;
}
