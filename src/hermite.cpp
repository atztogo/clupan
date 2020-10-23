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

    Class for obtaining Hermite normal form for derivative lattice

 ******************************************************************************/

#include "hermite.h"

Hermite::Hermite(const int& n_lattice, 
    const std::vector<imatrix>& rotate_matrix_array, const double& symprec){

    candidate_hermite_normal_matrix(n_lattice);
    unique_hnf(rotate_matrix_array, symprec);

}

Hermite::~Hermite(){}

void Hermite::candidate_hermite_normal_matrix(int n_lattice){

    // combination of diagonal terms
    // diagonal terms (a,c,f)
    std::vector<ivector> candidate_diagonal;
    for (int i = 1; i <= n_lattice; ++i){
        if ( n_lattice % i == 0) {
            int tmp1 = n_lattice / i;
            for (int j = 1; j <= tmp1; ++j){
                if ( tmp1 % j == 0) {
                    ivector list(3);
                    list(0) = i;
                    list(1) = j;
                    list(2) = tmp1/j;
                    candidate_diagonal.push_back(list);
                }
            }
        }
    }

    /************************  
      matrix 
      a 0 0
      b c 0 (0 <= b < c)
      d e f (0 <= d,e < f)  
     ************************/

    for (int i = 0; i < candidate_diagonal.size(); ++i){
        int a = candidate_diagonal[i](0);
        int c = candidate_diagonal[i](1);
        int f = candidate_diagonal[i](2);
        for (int b = 0; b < c; ++b){
            for (int d = 0; d < f; ++d){
                for (int e = 0; e < f; ++e){
                    imatrix candidate_matrix = ublas::zero_matrix<int>(3,3);
                    candidate_matrix(0,0) = a;
                    candidate_matrix(1,0) = b;
                    candidate_matrix(1,1) = c;
                    candidate_matrix(2,0) = d;
                    candidate_matrix(2,1) = e;
                    candidate_matrix(2,2) = f;
                    hermite_array.push_back(candidate_matrix);
                }
            }
        }
    }
}

void Hermite::unique_hnf(const std::vector<imatrix>& rotate_matrix_array,
        const double& symprec){

    std::set<int> remove;
    for (int i = 0; i < hermite_array.size(); ++i){
        dmatrix superlattice1;
        math::invert(hermite_array[i], superlattice1);
        for (int j = 0; j < rotate_matrix_array.size(); ++j){
            dmatrix prod_mat 
                = ublas::prod(superlattice1, rotate_matrix_array[j]);
            for (int k = 0; k < i; ++k){
                dmatrix prod_mat2 = ublas::prod(prod_mat, hermite_array[k]);
                int check = unimodular_check(prod_mat2, symprec);
                if (check == 1){
                    remove.insert(i);
                    goto loop;
                }
            }
        }
    loop:
        std::cout;
    }

    if (remove.size() > 0){
        std::set<int>::iterator it = remove.end();
        --it;
        while (it != remove.begin()){
            hermite_array.erase(hermite_array.begin() + *it);
            --it;
        }
        hermite_array.erase(hermite_array.begin() + *remove.begin());
    }

}

int Hermite::unimodular_check(const dmatrix& matrix, const double& symprec){

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

const std::vector<imatrix>& Hermite::get_hermite_normal_form_array() const{

    return hermite_array;

}
