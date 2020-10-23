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

    Class for labeling on derivative lattice

 ******************************************************************************/

#include "labeling.h"

Labeling::Labeling(){}

Labeling::~Labeling(){}

void Labeling::all_labeling(const std::vector<int>& n_atoms, 
    const std::vector<int>& n_type, 
    const std::vector<int>& n_atoms_derivative){

    std::vector<short> labeling_one_dim;

    // obtaining all configurations
    int add(0);
    for (int i = 0; i < n_type.size(); ++i){
        std::vector<short> labeling_single_lattice;
        Permutation perm;
        // in a whole range of compositions
        if (n_atoms_derivative.size() == 0){
            perm.permutation_label(n_atoms[i], n_type[i]);
            labeling_single_lattice = perm.get_permutation_label();
        }
        // in a fix composition
        else {
            std::vector<int> n_atoms_single;
            int sum(0);
            for (int j = 0; j < i; ++j){
                sum += n_type[i];
            }
            for (int j = 0; j < n_type[i]; ++j){
                n_atoms_single.push_back(n_atoms_derivative[sum+j]);
            }
            // n_atoms cannot work when n_atoms_single = 81 (Kuwabara)
            perm.permutation_label_fix_composition(n_atoms_single);
            labeling_single_lattice = perm.get_permutation_label();
        }
        for (int j = 0; j < labeling_single_lattice.size(); ++j){
            labeling_single_lattice[j] += add;
        }

        if (i > 0){
            int n_atom_labeling(0);
            for (int j = 0; j < i; ++j){
                n_atom_labeling += n_atoms[j];
            }
            combine_label(labeling_one_dim, labeling_single_lattice, 
                n_atom_labeling, n_atoms[i]);
        }
        else {
            labeling_one_dim = labeling_single_lattice;
        }
        add += n_type[i];
    }

    int n_all_atoms(0);
    for (int i = 0; i < n_type.size(); ++i){
        n_all_atoms += n_atoms[i];
    }

    // In a whole range of compositions, exchange duplicates are removed.

    if (n_atoms_derivative.size() == 0){
        n_type_check(labeling_one_dim, n_type, n_all_atoms);
    }
    else {
        int n(0);
        while (n < labeling_one_dim.size()){
            std::vector<short> tmp;
            for (int j = 0; j < n_all_atoms; ++j){
                tmp.push_back(labeling_one_dim[n]);
                ++n;
            }
            labeling.push_back(tmp);
        }
    }
}

void Labeling::combine_label(std::vector<short>& labeling_one_dim, 
    const std::vector<short>& labeling_single_lattice,
    const int& n_atom_labeling, const int& n_atom_labeling_single){

    std::vector<short> labeling_new;

    int n(0);
    while (n < labeling_one_dim.size()){
        std::vector<short> tmp;
        for (int j = 0; j < n_atom_labeling; ++j){
            tmp.push_back(labeling_one_dim[n]);
            ++n;
        }
        int n_s(0);
        while (n_s < labeling_single_lattice.size()){
            std::vector<short> tmp2(tmp);
            for (int j = 0; j < n_atom_labeling_single; ++j){
                tmp2.push_back(labeling_single_lattice[n_s]);
                ++n_s;
            }
            for (int j = 0; j < tmp2.size(); ++j){
                labeling_new.push_back(tmp2[j]);
            }
            tmp2.clear();
        }
    }
    labeling_one_dim = labeling_new;
    labeling_new.clear();

}

void Labeling::n_type_check(const std::vector<short>& label_candidate, 
    const std::vector<int>& n_type, const int& n_atom){

    labeling.clear();

    int n_type_sum(0);
    for (int i = 0; i < n_type.size(); ++i){
        n_type_sum += n_type[i];
    }

    int n_index(0);
    while (n_index < label_candidate.size()){
        std::vector<short> label;
        for (int i = 0; i < n_atom; ++i){
            label.push_back(label_candidate[n_index]);
            ++n_index;
        }
        std::vector<int> count(n_type_sum);
        for (int i = 0; i < n_atom; ++i){
            for (int j = 0; j < n_type_sum; ++j){
                if (label[i] == j){
                    count[j] += 1;
                }
            }
        }
        int n(0);
        for (int j = 0; j < n_type.size(); ++j){
            for (int k = 0; k < n_type[j] - 1; ++k){
                if (count[n] < count[n+1] 
                        or count[n] == 0 or count[n+1] == 0){
                    goto loop;
                }
                else if ( k == n_type[j] - 2 and j == n_type.size() - 1){
                    labeling.push_back(label);
                }
                ++n;
            }
            ++n;
        }
        loop:;
    }
}

void Labeling::permutation_check(
        const vector2_s& label_candidate, 
        const imatrix& product_table){

    std::set<std::vector<short> > unique_label_translation;

    for (int i = 0; i < label_candidate.size(); ++i){
        std::vector<short> label = label_candidate[i];
        std::set<std::vector<short> > duplicate_check;
        for (int j = 0; j < product_table.size1(); ++j){
            std::vector<short> trans;
            for (int k = 0; k < product_table.size2(); ++k){
                trans.push_back(label[product_table(j,k)]);
            }
            duplicate_check.insert(trans);
        }
        unique_label_translation.insert(*duplicate_check.begin());
    }

    labeling.clear();
    labeling.resize(unique_label_translation.size());
    std::copy(unique_label_translation.begin(), 
        unique_label_translation.end(), labeling.begin());

}

void Labeling::unimodular_sym_check(
    const imatrix& hnf, const imatrix& right,
    const std::vector<imatrix>& rotate_matrix_array,
    const std::vector<dvector>& trans_vector_array,
    std::vector<imatrix>& rotate_matrix_candidate,
    std::vector<dvector>& trans_vector_candidate,
    const double& symprec){

    dmatrix prod_mat = ublas::prod(hnf, right);
    dmatrix inverse_prod_mat;
    math::invert(prod_mat, inverse_prod_mat);

    for (int i = 0; i < rotate_matrix_array.size(); ++i){
        dmatrix prod_mat2 = ublas::prod(rotate_matrix_array[i], prod_mat);
        dmatrix prod_mat3 = ublas::prod(inverse_prod_mat, prod_mat2);
        int check = unimodular_check(prod_mat3, symprec);
        if (check == 1){
            rotate_matrix_candidate.push_back(rotate_matrix_array[i]);
            trans_vector_candidate.push_back(trans_vector_array[i]);
        }
    }
}

int Labeling::unimodular_check(
    const dmatrix& matrix, const double& symprec){

    int check(0);

    int count(0);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            double x = matrix(i,j);
            if (fabs(ceil(x) - x) < symprec or fabs(floor(x) - x) < symprec) {
                ++count;
            }
        }
    }
    if (count == 9){
        check = 1;
    }

    return check;

}

void Labeling::exchange_check(
    const vector2_s& label_original, 
    const std::vector<int>& n_type,
    const imatrix& permutation){

    vector2_i exchange = exchange_type(n_type);

    vector2_s labeling_output;
    for (int i = 0; i < label_original.size(); ++i){
        vector2_s label_new;
        for (int k = 0; k < exchange.size(); ++k){
            std::vector<short> tmp = label_original[i];
            for (int l = 0; l < exchange[k].size(); ++l){
                for (int m = 0; m < label_original[i].size(); ++m){
                    if (label_original[i][m] == l){
                        tmp[m] = exchange[k][l];
                    }
                }
            }
            label_new.push_back(tmp);
        }
        labeling.clear();
        permutation_check(label_new, permutation);
        labeling_output.push_back(labeling[0]);
    }

    labeling = labeling_output;
 
}

void Labeling::superlattice_check(
        const vector2_s& label_candidate, const dmatrix& product_table){

    std::set<std::vector<short> > unique_label_translation;
    for (int i = 0; i < label_candidate.size(); ++i){
        std::vector<short> label = label_candidate[i];
        std::set<std::vector<short> > duplicate_check;
        int n_order = product_table.size1();
        for (int j = 0; j < product_table.size1(); ++j){
            std::vector<short> trans;
            for (int k = 0; k < product_table.size2(); ++k){
                trans.push_back(label[product_table(j,k)]);
            }
            duplicate_check.insert(trans);
        }
        // superperiodic check
        if (duplicate_check.size() == n_order){
            unique_label_translation.insert(*duplicate_check.begin());
        }
    }

    labeling.clear();
    labeling.resize(unique_label_translation.size());
    std::copy(unique_label_translation.begin(), 
            unique_label_translation.end(), labeling.begin());

}

vector2_i Labeling::exchange_type(const std::vector<int>& n_type){

    vector2_i exchange;

    int add(0);
    for (int i = 0; i < n_type.size(); ++i){
        vector2_i exchange_single_lattice;
        Permutation perm;
        perm.permutation_without_repetition(n_type[i]);
        exchange_single_lattice = perm.get_permutation();
        for (int j = 0; j < exchange_single_lattice.size(); ++j){
            for (int k = 0; k < exchange_single_lattice[j].size(); ++k){
                exchange_single_lattice[j][k] += add;
            }
        }
        if (exchange.size() > 0){
            vector2_i exchange_new;
            for (int j = 0; j < exchange.size(); ++j){
                for (int k = 0; k < exchange_single_lattice.size(); ++k){
                    std::vector<int> tmp(exchange[j]);
                    tmp.insert(tmp.end(), exchange_single_lattice[k].begin(), 
                            exchange_single_lattice[k].end());
                    exchange_new.push_back(tmp);
                }
            }
            exchange = exchange_new;
        }
        else {
            exchange = exchange_single_lattice;
        }
        add += n_type[i];
    }

    return exchange;
}

const vector2_s& Labeling::get_label() const{
    return labeling;
}

