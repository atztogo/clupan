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

    Main program for finding special quasirandom structure (SQS) 
    using simulated annealing

******************************************************************************/

#include <iostream>
#include <vector>
#include <algorithm>

#include "input.h"
#include "output.h"
#include "mc_initialize.h"
#include "mc_binary.h"
#include "mc.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;

int main(){

    Input ip;
    ip.input_sqs();
    MC_initialize mc_init(ip);

    // checking whether system is binary or not
    std::vector<int> n_type = ip.get_n_type();
    int binary(0);
    if (n_type.size() == count(n_type.begin(), n_type.end(), 2)){
        binary = 1;
    }

    dvector correlation_array;
    std::vector<int> spin;
    if (binary == 1){
        MC_binary mc(ip, mc_init);
        mc.find_sqs(ip, mc_init);
        correlation_array = mc.get_correlation_array();
        spin = mc.get_spin();
    }
    else {
        MC mc(ip, mc_init);
        mc.find_sqs(ip, mc_init);
        correlation_array = mc.get_correlation_array();
        spin = mc.get_spin();
    }

    const dmatrix axis = mc_init.get_axis();
    const std::vector<int> n_atoms = mc_init.get_n_atoms();
    const std::vector<dvector> position = mc_init.get_position();
    std::vector<int> spin_array = ip.get_spin_array();
    const std::vector<int> cluster_index = mc_init.get_cluster_index();

    Output op;
    op.output_sqs(correlation_array, cluster_index,
            axis, n_atoms, position, n_type, spin, spin_array);

    return 0;
}
