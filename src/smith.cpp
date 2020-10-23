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

    Class for deriving Smith normal form

 ******************************************************************************/

#include "smith.h"

Smith::Smith(const imatrix& matrix){
 
    snf = matrix;
    row_matrix = ublas::identity_matrix<int>(3);
    column_matrix = ublas::identity_matrix<int>(3);

    if (snf.size1() != snf.size2()){
        std::cerr 
            << " error : only square matrix can be used for obtaining SNF. " 
            << std::endl;
        exit(8);
    }

    calc_smith_normal_form(snf.size1());
    calc_product_table();

}

Smith::~Smith(){}

void Smith::calc_smith_normal_form(int size){

    int start_index = snf.size1() - size;

    int min(1000);
    std::vector<int> min_index(2);
    for (int i = start_index; i < snf.size1(); ++i){
        for (int j = start_index; j < snf.size2(); ++j){
            if (abs(snf(i,j)) < min and snf(i,j) != 0){
                min_index[0] = i;
                min_index[1] = j;
                min = abs(snf(i,j));
            }
        }
    }

    imatrix permutation_row 
        = permutation_matrix(min_index[0], start_index, snf.size1());
    imatrix permutation_column 
        = permutation_matrix(min_index[1], start_index, snf.size1());

    snf = ublas::prod(permutation_row, snf);
    row_matrix = ublas::prod(permutation_row, row_matrix);

    snf = ublas::prod(snf, permutation_column);
    column_matrix = ublas::prod(column_matrix, permutation_column);

/*
    int check1(0), check2(0);
    for (int i = start_index + 1; i < snf.size1(); ++i){
        if (snf(i,start_index) % snf(start_index,start_index) != 0){
            ++check1;
        }
        if (snf(start_index,i) % snf(start_index,start_index) != 0){
            ++check2;
        }
    }
*/

    for (int i = start_index + 1; i < snf.size1(); ++i){
        int div = snf(i,start_index) / snf(start_index,start_index);
        imatrix add = addition_matrix(i, start_index, -div, snf.size1());
        snf = ublas::prod(add, snf);
        row_matrix = ublas::prod(add, row_matrix);
    }

    for (int i = start_index + 1; i < snf.size1(); ++i){
        int div = snf(start_index,i) / snf(start_index,start_index);
        imatrix add = addition_matrix(start_index, i, -div, snf.size1());
        snf = ublas::prod(snf, add);
        column_matrix = ublas::prod(column_matrix, add);
    }

    int zero_check = 0;
    for (int i = start_index + 1; i < snf.size1(); ++i){
        if (snf(i,start_index) != 0){
            zero_check = zero_check + 1;
            break;
        }
        if (snf(start_index, i) != 0){
            zero_check = zero_check + 1;
            break;
        }
    }

    if (zero_check == 0){
        int div_check(0);
        for (int i = start_index + 1; i < snf.size1(); ++i){
            for (int j = start_index + 1; j < snf.size1(); ++j){
                if (snf(i,j) % snf(start_index, start_index) != 0){
                    ++div_check;
                    imatrix add
                        = addition_matrix(start_index, i, 1, snf.size1());
                    snf = ublas::prod(add, snf);
                    row_matrix = ublas::prod(add, row_matrix);
                    break;
                }
            }
        }
        if (div_check == 0){
            size = size - 1;
        }
        if (size > 1){
            calc_smith_normal_form(size);
        }
        else {
            for (int i = 0; i < snf.size1(); ++i){
                if (snf(i,i) < 0){
                    imatrix multiply = multiply_matrix(i, -1, snf.size1());
                    snf = ublas::prod(multiply, snf);
                    row_matrix = ublas::prod(multiply, row_matrix);
                }
            }
        }
    }
    else {
        calc_smith_normal_form(size);
    }

}

void Smith::calc_product_table(){

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

    product_table = ublas::identity_matrix<int> (n_atom);
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
}



const imatrix& Smith::get_smith_normal_form() const{
    
    return snf;
    
}

const imatrix& Smith::get_row_matrix() const{
    
    return row_matrix;
    
}

const imatrix& Smith::get_column_matrix() const{
    
    return column_matrix;
    
}

const imatrix& Smith::get_product_table() const{
    
    return product_table;
    
}

imatrix Smith::permutation_matrix
    (const int& index1, const int& index2, const int& size){

    imatrix permutation = ublas::identity_matrix<int>(size);
    permutation(index1, index1) = 0;
    permutation(index2, index2) = 0;
    permutation(index1, index2) = 1;
    permutation(index2, index1) = 1;

    return permutation;
    
}

imatrix Smith::multiply_matrix
    (const int& index1, const int& multiplier, const int& size){

    imatrix multiply = ublas::identity_matrix<int>(size);
    multiply(index1, index1) = multiplier;

    return multiply;
    
}


imatrix Smith::addition_matrix
    (const int& index1, const int& index2, 
    const int& multiplier, const int& size){

    imatrix addition = ublas::identity_matrix<int>(size);
    addition(index1, index2) = multiplier;

    return addition;
    
}


