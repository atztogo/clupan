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

    Class for searching symmetry operations

 ******************************************************************************/

#include "sym.h"

Sym::Sym(const dmatrix& axis, const std::vector<dvector>& position,
        const std::vector<int>& type, const double symprec){

    get_sym(axis, position, type, symprec);

}

Sym::~Sym(){}

void Sym::get_sym(const dmatrix& axis, const std::vector<dvector>& position,
        const std::vector<int>& type, const double symprec){

    int natom(position.size());
    double lattice[3][3], position2[natom][3];
    int type2[natom];

    for (int i = 0; i < natom; ++i){
        for (int j = 0; j < 3; ++j){
            position2[i][j] = position[i](j);
        }
        type2[i] = type[i];
    }

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            lattice[i][j] = axis(i,j);
        }
    }
    
    int max_size = natom * 500;
    int rotation[max_size][3][3];
    double translation[max_size][3];

    int size;
    size = spg_get_symmetry(rotation, translation, max_size,
                        lattice, position2, type2, natom, symprec);
    
    for (int i = 0; i < size; ++i){
        imatrix tmpmat(3,3);
        dvector tmpvec(3);
        for (int j = 0; j < 3; ++j){
            tmpvec(j) = translation[i][j];
            for (int k = 0; k < 3; ++k){
                tmpmat(j,k) = rotation[i][j][k];
            }
        }
        rotate_matrix_array.push_back(tmpmat);
        trans_vector_array.push_back(tmpvec);
    }
//    for (int i = 0; i < rotate_matrix_array.size(); ++i){
 //       std::cout << rotate_matrix_array[i] << " "
  //          << trans_vector_array[i] << std::endl;
   // }

}
void Sym::lattice_based_primitive(const dmatrix& axis, 
    const std::vector<dvector>& position, 
    const std::vector<int>& type_lattice, const double& symprec,
    const std::vector<int> n_atoms, const dmatrix& axis_change, 
    const dmatrix& inverse_axis_change){

    std::vector<imatrix> rotate_matrix_array_primitive = rotate_matrix_array;
    std::vector<dvector> trans_vector_array_primitive = trans_vector_array;
    rotate_matrix_array.clear();
    trans_vector_array.clear();
/*
    // added on 2012/12/19
    imatrix identity_matrix = ublas::identity_matrix<int>(3);
    std::vector<dvector> trans_vector_remove;
    for (int i = 0; i < rotate_matrix_array_primitive.size(); ++i){
        imatrix rot = rotate_matrix_array_primitive[i];
        if (ublas::norm_inf(rot - identity_matrix) < symprec){
            trans_vector_remove.push_back(trans_vector_array_primitive[i]);
        }
    }

    std::vector<imatrix> rot_tmp;
    std::vector<dvector> trans_tmp;
    for (int i = 0; i < rotate_matrix_array_primitive.size(); ++i){
        for (int j = 1; j < trans_vector_remove.size(); ++j){
            dvector diff = trans_vector_array_primitive[i]
                - trans_vector_remove[j];
            if (ublas::inner_prod(diff,diff) < symprec){
                break;
            }
            else if (j == trans_vector_remove.size() - 1){
                rot_tmp.push_back(rotate_matrix_array_primitive[i]);
                trans_tmp.push_back(trans_vector_array_primitive[i]);
            }
        }
    }
    rotate_matrix_array_primitive = rot_tmp;
    trans_vector_array_primitive = trans_tmp;
    for (int i = 0; i < rotate_matrix_array_primitive.size(); ++i){
        std::cout << rotate_matrix_array_primitive[i] << " "
            << trans_vector_array_primitive[i] << std::endl;
    }
    ////////////////
*/   
    get_sym(axis, position, type_lattice, symprec);

    std::vector<dvector> pure_trans_vector;
    imatrix identity_matrix = ublas::identity_matrix<int>(3);
    for (int i = 0; i < rotate_matrix_array.size(); ++i){
        imatrix rot = rotate_matrix_array[i];
        if (ublas::norm_inf(rot - identity_matrix) < symprec){
            pure_trans_vector.push_back(trans_vector_array[i]);
        }
    }
    trans_vector_array.clear();

    for (int i = 0; i < rotate_matrix_array_primitive.size(); ++i){
        dvector trans_init = ublas::prod(inverse_axis_change,
                trans_vector_array_primitive[i]);
        dmatrix prodmat = ublas::prod(inverse_axis_change,
                rotate_matrix_array_primitive[i]);
        dmatrix prodmat2 = ublas::prod(prodmat, axis_change);
        for (int j = 0; j < pure_trans_vector.size(); ++j){
            double_rotate_matrix_array.push_back(prodmat2);
            trans_vector_array.push_back(trans_init + pure_trans_vector[j]);
        }
    }
//    std::cout << "double" << std::endl;
 //   for (int i = 0; i < double_rotate_matrix_array.size(); ++i){
  //      std::cout << double_rotate_matrix_array[i] << " "
   //         << trans_vector_array[i] << std::endl;
//    }
}

void Sym::supercell_symmetry(const dmatrix& axis_change){

    std::vector<imatrix> rotate_matrix_array_primitive = rotate_matrix_array;
    std::vector<dvector> trans_vector_array_primitive = trans_vector_array;
    rotate_matrix_array.clear();
    trans_vector_array.clear();

    dmatrix inverse_axis_change;
    math::invert(axis_change, inverse_axis_change);

    std::vector<dmatrix> double_rotate_matrix_array_primitive;
    for (int i = 0; i < rotate_matrix_array_primitive.size(); ++i){
        dmatrix prodmat = ublas::prod(inverse_axis_change,
                rotate_matrix_array_primitive[i]);
        dmatrix prodmat2 = ublas::prod(prodmat, axis_change);
        double_rotate_matrix_array_primitive.push_back(prodmat2);
    }


    for (int i = 0; i < double_rotate_matrix_array_primitive.size(); ++i){
        dvector trans_init = ublas::prod(inverse_axis_change,
                trans_vector_array_primitive[i]);
        for (int j = 0; j < round(axis_change(0,0)); ++j){
            for (int k = 0; k < round(axis_change(1,1)); ++k){
                for (int l = 0; l < round(axis_change(2,2)); ++l){
                    dvector trans(3);
                    trans(0) = j;
                    trans(1) = k;
                    trans(2) = l;
                    trans = ublas::prod(inverse_axis_change, trans);
                    dvector new_trans = trans_init + trans;
                    trans_vector_array.push_back(new_trans);
                    double_rotate_matrix_array.push_back
                        (double_rotate_matrix_array_primitive[i]);
                }
            }
        }
    }


}


const std::vector<imatrix>& Sym::get_rotate_matrix() const{

    return rotate_matrix_array;

}

const std::vector<dvector>& Sym::get_trans_vector() const{

    return trans_vector_array;

}

const std::vector<dmatrix>& Sym::get_double_rotate_matrix() const{

    return double_rotate_matrix_array;

}
