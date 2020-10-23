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

    Main program for calculating correlation functions

******************************************************************************/

#include <iostream>
#include <vector>
#include <set>

#include "input.h"
#include "output.h"
#include "error_check.h"
#include "sym.h"
#include "correlation.h"
#include "include/math.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef std::vector<std::vector<std::pair<ivector, int> > > vector2_pair;
typedef std::vector<std::vector<int> > vector2_i;

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

int main(int argc, char *argv[]){

    // reading input files
    Input ip;
    ip.input_correlation(argc, argv);

    const dmatrix axis_primitive(ip.get_axis_primitive());
    const std::vector<int> n_atoms_primitive(ip.get_n_atoms_primitive());
    const std::vector<dvector> position_primitive(ip.get_position_primitive());
    const std::vector<int> type_primitive(ip.get_type_primitive());

    const dmatrix axis(ip.get_axis());
    const std::vector<int> n_atoms(ip.get_n_atoms());
    const std::vector<dvector> position(ip.get_position());
    const std::vector<int> type(ip.get_type());

    const std::vector<int> n_type(ip.get_n_type());

    const dmatrix axis_change(ip.get_axis_change());
    const double symprec(ip.get_symprec());

    dmatrix inverse_axis_change;
    math::invert(axis_change, inverse_axis_change);

    // calculation of symmetric operations 
    // of lattice of structure (not structure)
    Sym symmetry(axis_primitive, position_primitive, type_primitive, symprec);
    symmetry.lattice_based_primitive(axis, position, type, 
        symprec, n_atoms, axis_change, inverse_axis_change);
    
    const std::vector<dmatrix> rotate_matrix_array
        (symmetry.get_double_rotate_matrix());
    const std::vector<dvector> trans_vector_array
        (symmetry.get_trans_vector());
    std::cout << "  number of symmetry operations = "
        << rotate_matrix_array.size() << std::endl;
    
    // unique clusters in primitive lattice 
    // -> permutations of clusters in "structure"
    const vector2_pair unique_cluster_array(ip.get_unique_cluster_array());
    Correlation correlation;
    std::vector<imatrix> permutation_array 
        = correlation.get_permutation_array(unique_cluster_array, 
        position_primitive, position, inverse_axis_change, 
        rotate_matrix_array, trans_vector_array, symprec);
/*
    for (int i = 0; i < permutation_array.size(); ++i){
        std::cout << " cluster " << i + 1<< std::endl;
        std::cout << permutation_array[i] << std::endl;
    }
*/
    // deriving cluster functions
    const std::vector<int> spin_array(ip.get_spin_array());

    std::vector<dmatrix> cluster_function_array;
    std::vector<int> cluster_number_array;
    if (n_type.size() == 1){
        cluster_function_array = correlation.get_cluster_function_array
            (spin_array, permutation_array, cluster_number_array);
    }
    else {
        vector2_i cluster_sublattice; 
        for (int i = 0; i < unique_cluster_array.size(); ++i){
            std::vector<int> sublattice;
            for (int j = 0; j < unique_cluster_array[i].size(); ++j){
                int atom = unique_cluster_array[i][j].second;
                for (int k = 0; k < n_atoms_primitive.size(); ++k){
                    atom -= n_atoms_primitive[k];
                    if (atom < 0){
                        sublattice.push_back(k);
                        break;
                    }
                }
            }
            cluster_sublattice.push_back(sublattice);
        }
/*        for (int i = 0; i < cluster_sublattice.size(); ++i){
            for (int j = 0; j < cluster_sublattice[i].size(); ++j){
                std::cout << cluster_sublattice[i][j];
            }
            std::cout << std::endl;
        }
*/
        cluster_function_array
            = correlation.get_cluster_function_array_multilattice
            (n_type, spin_array, cluster_sublattice, 
             permutation_array, cluster_number_array);
    }


    // calculation of correlation functions
    const std::vector<int> spin(ip.get_spin_structure());
    dvector correlation_array 
        = correlation.get_correlation_function(permutation_array, 
                cluster_function_array, spin);

    std::vector<int> n_clusters;
    for (int i = 0; i < permutation_array.size(); ++i){
        n_clusters.push_back(permutation_array[i].size1());
    }

    const char* file_name = ip.get_file_name();
    Output op;
    op.output_correlation(correlation_array, n_clusters, 
            cluster_function_array, cluster_number_array, file_name, n_type);

    return 0;

}

