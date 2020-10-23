/*****************************************************************************

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

    Class for obtaining permutation of atoms

 ******************************************************************************/

#include "permutation_atom.h"

Permutation_atom::Permutation_atom(){}

Permutation_atom::~Permutation_atom(){}

imatrix Permutation_atom::permutation_lattice_translation_multilattice(
    const imatrix& snf, const int& n_atom_primitive){

    imatrix permutation_lattice = permutation_lattice_translation(snf);

    int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);
    int n_atom = n_lattice * n_atom_primitive;

    imatrix permutation(permutation_lattice.size1(), n_atom);
    for (int j = 0; j < permutation_lattice.size1(); ++j){
        for (int k = 0; k < n_atom_primitive; ++k){
            for (int l = 0; l < n_lattice; ++l){
                int index = n_lattice * k + l;
                permutation(j,index)
                    = permutation_lattice(j, l) + n_lattice * k;
            }
        }
    }

    return permutation;

}

imatrix Permutation_atom::permutation_lattice_translation(const imatrix& snf){

    std::vector<imatrix> group_array;

    int n_atom(1);
    for (int i = 0; i < 3; ++i){
        imatrix group(snf(i,i), snf(i,i));
        n_atom *= snf(i,i);
        group_array.push_back(group);
    }

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < snf(i,i); ++j){
            for (int k = 0; k < snf(i,i); ++k){
                int index = j + k;
                if (index < snf(i,i)){
                    group_array[i](j,k) = index;
                }
                else {
                    group_array[i](j,k) = index - snf(i,i);
                }
            }
        }
    }

    imatrix product_table = ublas::identity_matrix<int> (n_atom);
    for (int i1 = 0; i1 < snf(0,0); ++i1){
        for (int i2 = 0; i2 < snf(0,0); ++i2){
            for (int j1 = 0; j1 < snf(1,1); ++j1){
                for (int j2 = 0; j2 < snf(1,1); ++j2){
                    for (int k1 = 0; k1 < snf(2,2); ++k1){
                        for (int k2 = 0; k2 < snf(2,2); ++k2){
                            int index1 = i1 * snf(1,1) * snf(2,2) 
                                + j1 * snf(2,2) + k1;
                            int index2 = i2 * snf(1,1) * snf(2,2) 
                                + j2 * snf(2,2) + k2;
                            product_table (index1,index2)
                                = group_array[0](i1,i2) * snf(1,1) * snf(2,2)
                                + group_array[1](j1,j2) * snf(2,2)
                                + group_array[2](k1,k2);
                        }
                    }
                }
            }
        }
    }

    return product_table;
}

imatrix Permutation_atom::permutation_rotation(
    const imatrix& snf, const imatrix& left, 
    const std::vector<imatrix>& rotate_matrix_array,
    const std::vector<dvector>& trans_vector_array,
    const std::vector<dvector>& position_primitive,
    const double& symprec){

    imatrix permutation_lattice 
        = permutation_lattice_translation_multilattice
        (snf, position_primitive.size());

    int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);
    int n_atom = n_lattice * position_primitive.size();

    imatrix position_lattice(3, n_lattice);
    for (int i = 0; i < n_lattice; ++i){
        ivector lattice(3);
        int atom_index;
        sequential_number_to_lattice_expression(i, snf, lattice, atom_index);
        for (int j = 0; j < 3; ++j){
            position_lattice(j, i) = lattice(j);
        }
    }

    dmatrix inverse_left;
    math::invert(left, inverse_left);

    std::vector<dvector> position_internal;
    for (int i = 0; i < position_primitive.size(); ++i){
        dvector pos = ublas::prod(left, position_primitive[i]);
        for (int j = 0; j < 3; ++j){
            while (pos(j) > 1 - symprec){
                pos(j) -= 1.0;
            }
            while (pos(j) < -symprec){
                pos(j) += 1.0;
            }
        }
        position_internal.push_back(pos);
    }
/*
    dmatrix axis(3,3);
    axis(0,0) = 4;
    axis(0,1) = -4;
    axis(0,2) = 4;
    axis(1,0) = -4;
    axis(1,1) = 0;
    axis(1,2) = 4;
    axis(2,0) = 0;
    axis(2,1) = -4;
    axis(2,2) = 0;
*/
    std::set<std::vector<int> > permutation_set;
    for (int i = 0; i < rotate_matrix_array.size(); ++i){
        dmatrix prod_mat = ublas::prod(left, rotate_matrix_array[i]);
        dmatrix prod_mat2 = ublas::prod(prod_mat, inverse_left);
        dmatrix prod_mat3 = ublas::prod(prod_mat2, position_lattice);
        dvector trans = ublas::prod(left, trans_vector_array[i]);

        std::vector<std::pair<dvector, dvector> > rotate_atoms;
        for (int atom = 0; atom < position_internal.size(); ++atom){
            dvector trans_pos 
                = ublas::prod(prod_mat2, position_internal[atom]) + trans;
            for (int j = 0; j < n_lattice; ++j){
                dvector lattice(3);
                for (int k = 0; k < 3; ++k){
                    lattice(k) =  prod_mat3(k,j);
                }
                std::pair<dvector, dvector> tmp_pair;
                tmp_pair.first = lattice;
                tmp_pair.second = trans_pos;
                rotate_atoms.push_back(tmp_pair);
            }
        }

        std::vector<int> permutation;
        for (int j = 0; j < rotate_atoms.size(); ++j){
            int permutation_index;
            dvector lattice = rotate_atoms[j].first;
            dvector trans_pos = rotate_atoms[j].second;
            ivector lattice_int(3);
            for (int k = 0; k < 3; ++k){
                while (trans_pos(k) > 1 - symprec){
                    trans_pos(k) -= 1.0;
                    lattice(k) += 1.0;
                }
                while (trans_pos(k) < -symprec){
                    trans_pos(k) += 1.0;
                    lattice(k) -= 1.0;
                }
                while (lattice(k) < -symprec){
                    lattice(k) += snf(k,k);
                }
                int lattice_integer = round(lattice(k));
                lattice_int(k) = lattice_integer % snf(k,k);
            }
            int atom_index;
            for (int k = 0; k < position_internal.size(); ++k){
                if (norm_inf(trans_pos - position_internal[k]) < symprec){
                    atom_index = k;
                    break;
                }
            }
            lattice_expression_to_sequential_number
                (lattice_int, atom_index, snf, permutation_index);
            permutation.push_back(permutation_index);
        }
        permutation_set.insert(permutation);
    }

    imatrix permutation_table(permutation_set.size() * n_lattice, n_atom);

    int n(0);
    std::set<std::vector<int> >::iterator it = permutation_set.begin();
    while (it != permutation_set.end()){
        std::vector<int> perm(*it);
        for (int j = 0; j < permutation_lattice.size1(); ++j){
            for (int k = 0; k < permutation_lattice.size2(); ++k){
                permutation_table(n,k) 
                    = perm[permutation_lattice(j,k)];
            }
            ++n;
        }
        ++it;
    }

    return permutation_table;

}

imatrix Permutation_atom::combine_table(const imatrix& rotation, 
    const imatrix& trans, const imatrix& snf){

    imatrix output(rotation.size1()*trans.size1(), trans.size2());

    int n(0);
    for (int i = 0; i < rotation.size1(); ++i){
        for (int k = 0; k < trans.size1(); ++k){
            for (int j = 0; j < rotation.size2(); ++j){
                ivector lattice_rotation, lattice_trans;
                int atom_index_rotation, atom_index_trans;
                sequential_number_to_lattice_expression
                    (rotation(i,j), snf, lattice_rotation, atom_index_rotation);
                sequential_number_to_lattice_expression
                    (trans(k,0), snf, lattice_trans, atom_index_trans);
                ivector new_lattice = lattice_trans + lattice_rotation;
                for (int l = 0; l < 3; ++l){
                    while (new_lattice(l) >= snf(l,l)){
                        new_lattice(l) -= snf(l,l);
                    }
                }
                int seq;
                lattice_expression_to_sequential_number
                    (new_lattice, atom_index_rotation, snf, seq);
                output(n, j) = seq;
            }
            ++n;
        }
    }

    return output;

}

imatrix Permutation_atom::cluster_rotation(
        const std::vector<dvector>& cluster_candidate,
        const std::vector<dmatrix>& rotate_matrix_array, 
        const std::vector<dvector>& trans_vector_array,
        const std::vector<dvector>& position,
        const double& symprec,
        const dmatrix& inverse_axis_change){

    dmatrix axis_change;
    math::invert(inverse_axis_change, axis_change);
/*
    dmatrix axis(3,3);
    axis(0,0) = 4;
    axis(0,1) = -4;
    axis(0,2) = 4;
    axis(1,0) = -4;
    axis(1,1) = 0;
    axis(1,2) = 4;
    axis(2,0) = 0;
    axis(2,1) = -4;
    axis(2,2) = 0;
*/
    std::vector<std::vector<dvector> > all_cluster;
    for (int i = 0; i < rotate_matrix_array.size(); ++i){
        std::vector<dvector> cluster;
        for (int j = 0; j < cluster_candidate.size(); ++j){
            dvector frac = cluster_candidate[j];
            dvector rotate_frac =
                ublas::prod(rotate_matrix_array[i], frac) 
                + trans_vector_array[i];
            cluster.push_back(rotate_frac);
        }
     
/////////////////////////////////////////////
        // NEEDED EXAMINATION
        // cluster translation to 

        std::vector<dvector> cluster_copy(cluster);
        dvector trans_cluster = ublas::zero_vector<double> (3);
        for (int k = 0; k < 3; ++k){
            while (cluster_copy[0](k) > 1 - symprec){
                cluster_copy[0](k) -= 1.0;
                trans_cluster(k) -= 1.0;
            }
            while (cluster_copy[0](k) < -symprec){
                cluster_copy[0](k) += 1.0;
                trans_cluster(k) += 1.0;
            }
        }
        if (cluster_copy.size() > 1){
            for (int k = 1; k < cluster_copy.size(); ++k){
                cluster_copy[k] += trans_cluster;
            }
        }

///////////////////////////////////////////

        if (all_cluster.size() == 0){
            all_cluster.push_back(cluster_copy);
        }
        else {
            for (int k = 0; k < all_cluster.size(); ++k){
                std::vector<dvector> cluster2 = all_cluster[k];
                for (int j = 0; j < cluster_copy.size(); ++j){
                    if (ublas::norm_inf(cluster_copy[j] - cluster2[j]) 
                        > symprec){
                        break;
                    }
                    else if (j == cluster_copy.size() - 1){
                        goto loop;
                    }
                }
                if (k == all_cluster.size() - 1){
                    all_cluster.push_back(cluster_copy);
                }
            }
            loop:;
        }
    }

    imatrix permutation(all_cluster.size(), cluster_candidate.size());

    for (int i = 0; i < all_cluster.size(); ++i){
        for (int j = 0; j < all_cluster[i].size(); ++j){
            dvector pos = all_cluster[i][j];
            for (int k = 0; k < pos.size(); ++k){
                while (pos(k) > 1 - symprec){
                    pos(k) -= 1.0;
                }
                while (pos(k) < -symprec){
                    pos(k) += 1.0;
                }
            }
            for (int k = 0; k < position.size(); ++k){
                if (ublas::norm_inf(pos - position[k]) < symprec){
                    permutation(i, j) = k;
                    break;
                }
                else if (k == position.size() - 1){
                    std::cerr << " error: position of rotated cluster ";
                    std::cerr << "cannot find in structure" << std::endl;
                    std::cerr << " index of rotational operation" << std::endl;
                    std::cerr << rotate_matrix_array[i] << std::endl;
                    std::cerr << trans_vector_array[i] << std::endl;
                    std::cerr << " original cluster" << std::endl;
                    for (int l = 0; l < cluster_candidate.size(); ++l){
                        std::cerr << cluster_candidate[l] << " ";
                    }
                    std::cerr << std::endl;
                    std::cerr << " rotated cluster" << std::endl;
                    for (int l = 0; l < all_cluster[i].size(); ++l){
                        std::cerr << all_cluster[i][l] << " ";
                    }
                    std::cerr << std::endl;
                    exit(8);
                }
            }
        }
    }
    return permutation;

}
