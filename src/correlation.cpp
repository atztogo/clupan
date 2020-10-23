/***************************************************************************

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
        
        Class of functions for calculating correlation functions

***************************************************************************/

#include "correlation.h"

Correlation::Correlation(){}
Correlation::~Correlation(){}

std::vector<imatrix> Correlation::get_permutation_array(
    const vector2_pair& unique_cluster_array,
    const std::vector<dvector> position_primitive,
    const std::vector<dvector> position,
    const dmatrix& inverse_axis_change,
    const std::vector<dmatrix>& rotate_matrix_array,
    const std::vector<dvector>& trans_vector_array,
    const double& symprec){

    // fractional coordinate of unique clusters in the basis of "structure"
    std::vector<std::vector<dvector> > unique_cluster_array_in_structure;
    unique_cluster_array_in_structure.resize(unique_cluster_array.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < unique_cluster_array.size(); ++i){
        std::vector<dvector> cluster_in_structure;
        for (int j = 0; j < unique_cluster_array[i].size(); ++j){
            ivector lattice = unique_cluster_array[i][j].first;
            dvector pos = position_primitive[unique_cluster_array[i][j].second];
            dvector frac = ublas::prod(inverse_axis_change, lattice + pos);
            cluster_in_structure.push_back(frac);
        }
        //        unique_cluster_array_in_structure.push_back(cluster_in_structure);
        unique_cluster_array_in_structure[i] = cluster_in_structure;
    }

    // all equivalent clusters rotated by symmetric operations
    std::vector<imatrix> permutation_array;
    permutation_array.resize(unique_cluster_array_in_structure.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < unique_cluster_array_in_structure.size(); ++i){
        Permutation_atom perm_atom;
        imatrix permutation = perm_atom.cluster_rotation
            (unique_cluster_array_in_structure[i], rotate_matrix_array,
             trans_vector_array, position, symprec, inverse_axis_change);
        //        permutation_array.push_back(permutation);
        permutation_array[i] = permutation;
    }

    return permutation_array;

}

std::vector<dmatrix> Correlation::get_cluster_function_array(
    const std::vector<int>& spin_array,
    std::vector<imatrix>& permutation_array,
    std::vector<int>& cluster_number_array){

    Hamiltonian hamiltonian(spin_array);
    dmatrix onsite_function = hamiltonian.get_onsite_function();

    std::vector<imatrix> permutation_array_output;
    std::vector<dmatrix> cluster_function_array;

    for (int i = 0; i < permutation_array.size(); ++i){
        int n_atom_cluster = permutation_array[i].size2();
        int n_type = spin_array.size();
        Permutation perm;
        perm.permutation_with_repetition(n_atom_cluster, n_type -1);
        vector2_i permutation_cluster_function = perm.get_permutation();

        for (int j = 0; j < permutation_cluster_function.size(); ++j){
            dmatrix function(n_atom_cluster, n_type);
            for (int k = 0; k < n_atom_cluster; ++k){
                for (int l = 0; l < n_type; ++l){
                    function(k,l) = onsite_function
                        (permutation_cluster_function[j][k] + 1, l);
                }
            }
            cluster_function_array.push_back(function);
            permutation_array_output.push_back(permutation_array[i]);
            cluster_number_array.push_back(i);
        }
    }

    permutation_array = permutation_array_output;

    return cluster_function_array;

}

std::vector<dmatrix> Correlation::get_cluster_function_array_multilattice(
    const std::vector<int>& n_type,
    const std::vector<int>& spin_array,
    const vector2_i& cluster_sublattice,
    std::vector<imatrix>& permutation_array,
    std::vector<int>& cluster_number_array){

    // onsite functions for all sublattices
    std::vector<dmatrix> onsite_function_array;
    int n(0);
    for (int i = 0; i < n_type.size(); ++i){
        std::vector<int> spin_array_single;
        for (int j = 0; j < n_type[i]; ++j){
            spin_array_single.push_back(spin_array[n]);
            ++n;
        }
        Hamiltonian hamiltonian(spin_array_single);
        dmatrix onsite_function = hamiltonian.get_onsite_function();
        onsite_function_array.push_back(onsite_function);
    }

    std::vector<imatrix> permutation_array_output;
    std::vector<dmatrix> cluster_function_array;

    int n_type_sum(0);
    for (int i = 0; i < n_type.size(); ++i){
        n_type_sum += n_type[i];
    }
    int max_type = *max_element(n_type.begin(), n_type.end());

    for (int i = 0; i < cluster_sublattice.size(); ++i){
        int n_atom_cluster = cluster_sublattice[i].size();
        std::vector<int> n_type_array;
        for (int j = 0; j < cluster_sublattice[i].size(); ++j){
            n_type_array.push_back(n_type[cluster_sublattice[i][j]] - 1);
        }
        Permutation perm;
        perm.permutation_with_repetition_multilattice(n_type_array);
        vector2_i permutation_cluster_function = perm.get_permutation();
        for (int j = 0; j < permutation_cluster_function.size(); ++j){
            dmatrix function 
                = ublas::zero_matrix<double> (n_atom_cluster, max_type);
            for (int k = 0; k < n_atom_cluster; ++k){
                for (int l = 0; l < n_type[cluster_sublattice[i][k]]; ++l){
                    function(k,l) 
                    = onsite_function_array[cluster_sublattice[i][k]]
                        (permutation_cluster_function[j][k] + 1, l);
                }
            }
            cluster_function_array.push_back(function);
            permutation_array_output.push_back(permutation_array[i]);
            cluster_number_array.push_back(i);
        }
    }

    permutation_array = permutation_array_output;

    return cluster_function_array;

}


dvector Correlation::get_correlation_function(
    const std::vector<imatrix>& permutation_array,
    const std::vector<dmatrix>& cluster_function_array,
    const std::vector<int>& spin){

    dvector correlation_array(permutation_array.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < permutation_array.size(); ++i){
        dmatrix cluster_function = cluster_function_array[i];
        double correlation(0.0);
        for (int j = 0; j < permutation_array[i].size1(); ++j){
            std::vector<int> spin_values;
            for (int k = 0; k < permutation_array[i].size2(); ++k){
                spin_values.push_back(spin[permutation_array[i](j,k)]);
            }
            double product = calc_product(spin_values, cluster_function);
            correlation += product / double(permutation_array[i].size1());
        }
        correlation_array(i) = correlation;
    }

    return correlation_array;
}

double Correlation::calc_product
(std::vector<int>& spin_values, dmatrix& function){

    double product(1.0);
    for (int i = 0; i < function.size1(); ++i){
        double sum(0.0);
        for (int j = 0; j < function.size2(); ++j){
            double tmp_prod(1.0);
            for (int k = 0; k < j; ++k){
                tmp_prod *= spin_values[i];
            }
            sum += function(i,j) * tmp_prod;
        }
        product *= sum;
    }

    return product;
}

