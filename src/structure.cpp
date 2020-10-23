/******************************************************************************

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

	    Class for structure
	
*****************************************************************************/

#include "structure.h"

Structure::Structure
(const dmatrix& axis_input, const std::vector<dvector>& position_input,
 const std::vector<int>& n_atoms_input)
:axis(axis_input), position(position_input), n_atoms(n_atoms_input){}

Structure::~Structure(){}

double Structure::get_volume(){

    ublas::matrix_column<dmatrix> axis1(axis, 0);
    ublas::matrix_column<dmatrix> axis2(axis, 1);
    ublas::matrix_column<dmatrix> axis3(axis, 2);
    volume = ublas::inner_prod(axis1, math::cross_prod(axis2, axis3));

    return volume;
}

dmatrix Structure::get_reciprocal_axis(){

    get_volume();

    const double pi(3.14159265358979323846);
    ublas::matrix_column<dmatrix> axis1(axis, 0);
    ublas::matrix_column<dmatrix> axis2(axis, 1);
    ublas::matrix_column<dmatrix> axis3(axis, 2);
    const double coeff_reciprocal(2 * pi / volume);
    dvector rec_axis1 = coeff_reciprocal * math::cross_prod(axis2,axis3);
    dvector rec_axis2 = coeff_reciprocal * math::cross_prod(axis3,axis1);
    dvector rec_axis3 = coeff_reciprocal * math::cross_prod(axis1,axis2);

    dmatrix reciprocal_axis(3,3);
    for (int i = 0; i < 3; ++i){
        reciprocal_axis(i,0) = rec_axis1(i);
        reciprocal_axis(i,1) = rec_axis2(i);
        reciprocal_axis(i,2) = rec_axis3(i);
    }
    return reciprocal_axis;
}

int Structure::supercell
(const imatrix& supercell_matrix){

    // axes of supercell
    axis = ublas::prod(axis, supercell_matrix);

    // fractional coordinates in supercell
    imatrix snf, left;
    dmatrix inverse_supercell_matrix, inverse_left;

    //imatrix supercell_imatrix;
    //eigen_matrix_to_boost_matrix(supercell_matrix, supercell_imatrix);
    Smith smith(supercell_matrix);
    snf = smith.get_smith_normal_form();
    left = smith.get_row_matrix();
    math::invert(supercell_matrix, inverse_supercell_matrix);
    math::invert(left, inverse_left);

    // SNF = left * H * right
    // cart = A * x 
    //      = A * H * right * z  
    // --> left * A^(-1) * cart = left * x = S*z

    // S*z = labels of lattice points
    //  For example S*z = [0,0], [0,1], [1,0], [1,1] 
    //  when S = matrix((2,0), (0,2))

    // S*z of lattice points --> x (frac. on basis A)
    // left * x = S * z (frac. coord. on AH)
    // x = left^(-1) * S * z

    int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);
    imatrix position_lattice(3, n_lattice);
    for (int i = 0; i < n_lattice; ++i){
        ivector lattice; int atom_index;
        sequential_number_to_lattice_expression(i, snf, lattice, atom_index);
        for (int j = 0; j < 3; ++j){
            position_lattice(j, i) = lattice(j);
        }
    }
    dmatrix origin_frac = ublas::prod(inverse_left, position_lattice);

    // y = internal positions
    // right * z = H^(-1) * (x+y)
    std::vector<dvector> position_new;
    dvector pos_tmp, new_pos;
    for (int i = 0; i < position.size(); ++i){
        for (int j = 0; j < n_lattice; ++j){
            pos_tmp = ublas::column(origin_frac, j) + position[i];
            new_pos = ublas::prod(inverse_supercell_matrix, pos_tmp);
            for (int k = 0; k < 3; ++k){
                new_pos(k) -= floor(new_pos(k));
            }
            position_new.push_back(new_pos);
        }
    }
    position = position_new;

    // number of atoms in supercell
    for (int i = 0; i < n_atoms.size(); ++i){
        n_atoms[i] *= n_lattice;
    }
}

const dmatrix& Structure::get_axis() const{
    return axis;
}
const std::vector<int>& Structure::get_n_atoms() const{
    return n_atoms;
}
const std::vector<dvector>& Structure::get_position() const{
    return position;
}
