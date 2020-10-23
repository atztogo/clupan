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

    Class for obtaining permutation

 ******************************************************************************/

#include "permutation.h"

Permutation::Permutation(){}

Permutation::~Permutation(){}

void Permutation::permutation_label
    (const int& n_atom, const int& n_type){

    long n_comb = pow (n_type, n_atom);
    for (long i = 0; i < n_comb; ++i){
        short comb[n_atom];
        for (int j = 0; j < n_atom; ++j){
            comb[j] = 0;
        }
        int digit(n_atom-1);
        long n10 = i;
        while (n10 != 0){
            short a = n10 % n_type;
            comb[digit] = a;
            n10 /= n_type;
            --digit;
        }
        for (int j = 0; j < n_atom; ++j){
            permutation_one_dimension.push_back(comb[j]);
        }
    }
}

void Permutation::permutation_label_fix_composition
    (const std::vector<int>& n_atoms_derivative){

    int n_atom(0);
    for (int i = 0; i < n_atoms_derivative.size(); ++i){
        n_atom += n_atoms_derivative[i];
    }
    int n_type = n_atoms_derivative.size();
    
    // n_comb cannot work when n=2^81 (Kuwabara)
    long n_comb = pow (n_type, n_atom);
    for (long i = 0; i < n_comb; ++i){
        short comb[n_atom];
        for (int j = 0; j < n_atom; ++j){
            comb[j] = 0;
        }
        int digit(n_atom-1);
        long n10 = i;
        while (n10 != 0){
            short a = n10 % n_type;
            comb[digit] = a;
            n10 /= n_type;
            --digit;
        }
        int check = composition_check(comb, n_atoms_derivative);
        if (check == 1){
            for (int j = 0; j < n_atom; ++j){
                permutation_one_dimension.push_back(comb[j]);
            }
        }
    }
}

int Permutation::composition_check(short comb[ ], 
    const std::vector<int>& n_atoms_derivative){

    int check(0);

    int n_atom(0);
    for (int i = 0; i < n_atoms_derivative.size(); ++i){
        n_atom += n_atoms_derivative[i];
    }

    std::vector<int> count;
    for (int i = 0; i < n_atoms_derivative.size(); ++i){
        count.push_back(0);
    }
    for (int i = 0; i < n_atom; ++i){
        count[comb[i]] += 1;
    }
    for (int i = 0; i < n_atoms_derivative.size(); ++i){
        if (n_atoms_derivative[i] != count[i]){
            break;
        }
        else if (i == n_atoms_derivative.size() - 1){
            check = 1;
        }
    }

    return check;

}

void Permutation::permutation_without_repetition(const int& n_type){

    gsl_permutation * p = gsl_permutation_alloc (n_type);
    gsl_permutation_init (p);

    do {
        std::vector<int> element;
        for (int i = 0; i < n_type; ++i){
            element.push_back(gsl_permutation_get(p, i));
        }
        permutation.push_back(element);
    }
    while (gsl_permutation_next(p) == GSL_SUCCESS);
    gsl_permutation_free (p);

}

void Permutation::permutation_with_repetition
    (const int& n_atom, const int& n_type){

    long n_comb = pow (n_type, n_atom);
    std::vector<int> zero;
    for (int i = 0; i < n_atom; ++i){
        zero.push_back(0);
    }

    for (long i = 0; i < n_comb; ++i){
        std::vector<int> comb(zero);
        int digit(n_atom-1);
        long n10 = i;
        while (n10 != 0){
            int a = n10 % n_type;
            comb[digit] = a;
            n10 /= n_type;
            --digit;
        }
        permutation.push_back(comb);
    }

}

void Permutation::permutation_with_repetition_multilattice(
    const std::vector<int>& n_type_array){

    int n_comb(1);
    for (int i = 0; i < n_type_array.size(); ++i){
        n_comb *= n_type_array[i];
    }

    for (int i = 0; i < n_comb; ++i){
        std::vector<int> comb(n_type_array.size());
        int quotient(i);
        for (int j = 0; j < n_type_array.size(); ++j){
            int index = n_type_array.size()-j-1;
            comb[index] = quotient % n_type_array[index];
            quotient /= n_type_array[index];
        }
        permutation.push_back(comb);
    }
}

const std::vector<short>& Permutation::get_permutation_label() const{
    return permutation_one_dimension;
}
const vector2_i& Permutation::get_permutation() const{
    return permutation;
}
