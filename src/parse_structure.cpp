/******************************************************************************

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

	    Class for parsing POSCAR
	
*****************************************************************************/

#include "parse_structure.h"

Parse_structure::Parse_structure(const char *filename){

    Parse_input input(filename);

    dvector axis1, axis2, axis3;

    input.assignNeed("axis1", axis1);
    input.assignNeed("axis2", axis2);
    input.assignNeed("axis3", axis3);

    dmatrix axis_tmp(3,3);
    for (int i = 0; i < 3; ++i){
        axis_tmp(i,0) = axis1(i);
        axis_tmp(i,1) = axis2(i);
        axis_tmp(i,2) = axis3(i);
    }
    axis = axis_tmp;

    input.assignNeed("n_atoms", n_atoms);

    int count(1);
    for (int i = 0; i < n_atoms.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            dvector pos;
            std::string str("position"), str2;
            std::stringstream ss;
            ss << count;
            ss >> str2;
            str += str2;
            input.assignNeed(str.c_str(), pos);
            for (int k = 0; k < 3; ++k){
                pos(k) = modify_position(pos(k));
            }
            position.push_back(pos);
            ++count;
        }
    }
    
    for (int i = 0; i < n_atoms.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            type.push_back(i);
        }
    }
}

Parse_structure::~Parse_structure(){}

const dmatrix& Parse_structure::get_axis() const{

    return axis;
}

const std::vector<int>& Parse_structure::get_num_atoms() const{

    return n_atoms;
}

const std::vector<dvector>& Parse_structure::get_position() const{

    return position;
}

const std::vector<int>& Parse_structure::get_type() const{

    return type;
}

double Parse_structure::modify_position(double x) const{

    double x_mod;

    if (x > 1 - 1e-10)
        x_mod = x - floor(x);
    else if (x < - 1e-10)
        x_mod = x - floor(x);
    else
        x_mod = x;

    return x_mod;

}

double Parse_structure::get_volume(){


    ublas::matrix_column<dmatrix> axis1(axis, 0);
    ublas::matrix_column<dmatrix> axis2(axis, 1);
    ublas::matrix_column<dmatrix> axis3(axis, 2);

    volume = ublas::inner_prod(axis1, math::cross_prod(axis2, axis3));

    return volume;

}

dmatrix Parse_structure::get_reciprocal_axis(){

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

