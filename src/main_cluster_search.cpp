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

    Main program for searching for unique clusters 

******************************************************************************/

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "input.h"
#include "output.h"
#include "sym.h"
#include "combination.h"
#include "permutation_atom.h"
#include "include/label.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef std::vector<std::vector<int> > vector2_i;

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<int> imatrix;
typedef ublas::vector<int> ivector;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

int main(){

    // reading input files
    Input ip;
    ip.input_cluster();
	int max_cluster = ip.get_max_cluster();
    std::vector<double> trunc_distance = ip.get_trunc_distance();
    double symprec = ip.get_symprec();
    int n_sublattice = ip.get_n_sublattice();

    dmatrix axis_primitive = ip.get_axis_primitive();
    std::vector<int> n_atoms = ip.get_n_atoms_primitive();
    std::vector<dvector> position_primitive = ip.get_position_primitive();
    std::vector<int> type = ip.get_type_primitive();

    dvector axis_primitive_length(3);
    for (int i = 0; i < 3; ++i){
        ublas::matrix_column<dmatrix> axis_component(axis_primitive, i);
        axis_primitive_length(i) 
            = sqrt(ublas::inner_prod(axis_component, axis_component));
    }

    // calculation of symmetric operations
    Sym symmetry(axis_primitive, position_primitive, type, symprec);
    const std::vector<imatrix> 
        rotate_matrix_array(symmetry.get_rotate_matrix());
    const std::vector<dvector> 
        trans_vector_array(symmetry.get_trans_vector());
    std::cout << "  number of symmetry operations = " 
        << rotate_matrix_array.size() << std::endl;

    // positions on isub
    std::vector<dvector> position_simulation =  ip.get_position_simulation();
    int n_atom_primitive = position_simulation.size();

    // calculation of cell size for identifying unique clusters

    imatrix cell_expand = ublas::identity_matrix<int>(3);
    double max_trunc_distance = 
        *std::max_element(trunc_distance.begin(), trunc_distance.end());
    for (int i = 0; i < 3; ++i){
        double dis(axis_primitive_length(i));
        cell_expand(i,i) = int(3 * max_trunc_distance / dis) + 1;
    }
//    for (int i = 0; i < 3; ++i){
 //       cell_expand(i,i) = 20;
  //  }

    int n_lattice = cell_expand(0,0) * cell_expand(1,1) * cell_expand(2,2);
    int n_atom = n_lattice * n_atom_primitive;

    // obtaining permutation of atoms by symmetric operations
    Permutation_atom perm_atom;
    imatrix permutation_table = perm_atom.permutation_rotation
        (cell_expand, ublas::identity_matrix<int>(3), 
         rotate_matrix_array, trans_vector_array,
         position_simulation, symprec);

    // calculation of distance between two atoms
    ublas::symmetric_matrix<double> distance_array(n_atom, n_atom);
    int atom1(0);
    for (int i = 0; i < n_atom_primitive; ++i){
        for (int j = 0; j < n_lattice; ++j){
            for (int atom2 = 0; atom2 < atom1; ++atom2){
                ivector lattice1(3), lattice2(3);
                int index1, index2;
                sequential_number_to_lattice_expression
                    (atom1, cell_expand, lattice1, index1);
                sequential_number_to_lattice_expression
                    (atom2, cell_expand, lattice2, index2);
                dvector diff_pos = lattice2 + position_simulation[index2] 
                    - (lattice1 + position_simulation[index1]);
                for (int k = 0; k < 3; ++k){
                    if (diff_pos(k) > double(cell_expand(k,k)) / 2.0){
                        diff_pos(k) -= cell_expand(k,k); 
                    }
                    else if (diff_pos(k) < -double(cell_expand(k,k)) / 2.0){
                        diff_pos(k) += cell_expand(k,k);
                    }
                }
                dvector diff_pos_cart = ublas::prod(axis_primitive, diff_pos);
                distance_array(atom1,atom2) 
                    = sqrt(ublas::inner_prod(diff_pos_cart, diff_pos_cart));
            }
            ++atom1;
        }

    }
/*
    for (int i = 0; i < permutation_table.size1(); ++i){
        int atom1 = permutation_table(i,2000);
        int atom2 = permutation_table(i,3113);
        if (atom1 == 2000 or atom2 == 2000 ){
        std::cout << atom1 << " " << atom2 << " " 
            << distance_array(atom1, atom2) << std::endl;

        }
    }
*/
    // obtaining unique clusters
    std::cout << "  number of empty cluster = 1" << std::endl;
    vector2_i old_unique_cluster, all_unique_cluster;
    for (int n_cluster = 1; n_cluster < max_cluster + 1; ++n_cluster){
        double trunc_distance_tmp(0.0);
        if (n_cluster == 2){
            trunc_distance_tmp = trunc_distance[0];
        }
        else if (n_cluster == 3){
            trunc_distance_tmp = trunc_distance[1];
        }
        else if (n_cluster == 4){
            trunc_distance_tmp = trunc_distance[2];
        }
        else if (n_cluster > 4){
            trunc_distance_tmp = trunc_distance[3];
        }

        vector2_i combination_array;
        if (n_cluster == 1){
            for (int i = 0; i < n_atom; ++i){
                std::vector<int> tmp;
                tmp.push_back(i);
                combination_array.push_back(tmp);
            }
        }
        else if (n_cluster > 1){
            for (int i = 0; i < old_unique_cluster.size(); ++i){
                Combination comb;

                if (n_cluster > 2){
                    imatrix atom_comb 
                        = comb.get_combination(old_unique_cluster[i].size(), 2);
                    for (int j = 0; j < atom_comb.size1(); ++j){
                        int atom1 = old_unique_cluster[i][atom_comb(j,0)];
                        int atom2 = old_unique_cluster[i][atom_comb(j,1)];
                        double dis = distance_array(atom1,atom2);
                        if (dis > trunc_distance_tmp){
                            goto loop;
                        }
                    }
                }

                for (int j = 0; j < n_atom; ++j){
                    std::vector<int> tmp(old_unique_cluster[i]);
                    if (find(old_unique_cluster[i].begin(), 
                        old_unique_cluster[i].end(), j) 
                        == old_unique_cluster[i].end()){
                        for (int k = 0; k < old_unique_cluster[0].size(); ++k){
                            double dis = distance_array(j, tmp[k]);
                            if (dis > trunc_distance_tmp){
                                break;
                            }
                            else if (k == old_unique_cluster[0].size() - 1){
                                tmp.push_back(j);
                                combination_array.push_back(tmp);
                            }
                        }
                    }
                }
                loop:;
            }
        }

        std::set<std::vector<int> > unique_cluster_set;
        for (int i = 0; i < combination_array.size(); ++i){
            std::set<std::vector<int> > duplicate_check;
            for (int j = 0; j < permutation_table.size1(); ++j){
                std::vector<int> permutation_combination;
                for (int k = 0; k < combination_array[i].size(); ++k){
                    permutation_combination.push_back
                        (permutation_table(j, combination_array[i][k]));
                }
                std::sort(permutation_combination.begin(), 
                    permutation_combination.end());
                duplicate_check.insert(permutation_combination);
            }
            unique_cluster_set.insert(*duplicate_check.begin());
        }

        old_unique_cluster.clear();
        old_unique_cluster.resize(unique_cluster_set.size());
        std::copy(unique_cluster_set.begin(), unique_cluster_set.end(), 
            old_unique_cluster.begin());
        std::cout << "  number of " << n_cluster << "-body clusters = " 
            << old_unique_cluster.size() << std::endl;

        if (n_cluster == 2){
            std::multimap<double, std::vector<int> > sort;
            for (int i = 0; i < old_unique_cluster.size(); ++i){
                std::pair<double, std::vector<int> > tmp_pair;
                int atom1 = old_unique_cluster[i][0];
                int atom2 = old_unique_cluster[i][1];
                tmp_pair.first = distance_array(atom1,atom2);
                tmp_pair.second = old_unique_cluster[i];
                sort.insert(tmp_pair);
            }
            std::multimap<double, std::vector<int> >
                ::iterator it = sort.begin();
            while (it != sort.end()){
                all_unique_cluster.push_back((*it).second);
                ++it;
            }
        }
        else {
            for (int i = 0; i < old_unique_cluster.size(); ++i){
                all_unique_cluster.push_back(old_unique_cluster[i]);
            }
        }
    }

    // correction of clusters by periodic boundary condition
    std::vector<std::vector<std::pair<ivector, int> > > 
        all_unique_cluster_output;
    for (int i = 0; i < all_unique_cluster.size(); ++i){
        std::vector<std::pair<ivector, int> > cluster;
        for (int j = 0; j < all_unique_cluster[i].size(); ++j){
            std::pair<ivector, int> atom;
            sequential_number_to_lattice_expression
                (all_unique_cluster[i][j], cell_expand, 
                atom.first, atom.second);
            cluster.push_back(atom);
        }
        ivector lattice1 = cluster[0].first;
        int index1 = cluster[0].second;
        for (int j = 1; j < cluster.size(); ++j){
            ivector lattice2 = cluster[j].first;
            int index2 = cluster[j].second;
            dvector diff_pos = lattice2 + position_simulation[index2] 
                - (lattice1 + position_simulation[index1]);
            for (int k = 0; k < 3; ++k){
                if (diff_pos(k) > double(cell_expand(k,k)) / 2.0){
                    cluster[j].first(k) -= cell_expand(k,k); 
                }
                else if (diff_pos(k) < -double(cell_expand(k,k)) / 2.0){
                    cluster[j].first(k) += cell_expand(k,k);
                }
            }
        }
        all_unique_cluster_output.push_back(cluster);
    }

    std::cout << "  number of all clusters = " 
        << all_unique_cluster_output.size() + 1 << std::endl;

    // output
    Output op;
    op.output_cluster_out(all_unique_cluster_output, all_unique_cluster, 
        axis_primitive, n_atoms, position_primitive, 
        n_sublattice, distance_array);

    return 0;

}
