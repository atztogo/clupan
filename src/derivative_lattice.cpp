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

    Class for derivative lattice

 ******************************************************************************/

#include "derivative_lattice.h"

Derivative_lattice::Derivative_lattice(
    const imatrix& hnf_input, const imatrix& snf_input,
    const imatrix& left_input, const imatrix& right_input,
    const imatrix& product_table_input, const int& snf_index_input, 
    const vector2_s& label_input)
    :hnf(hnf_input), snf(snf_input), left(left_input), right(right_input), 
    product_table(product_table_input), snf_index(snf_index_input), 
    label(label_input){

}

Derivative_lattice::~Derivative_lattice(){}

const imatrix& Derivative_lattice::get_hermite_normal_form() const{
    
    return hnf;
    
}

const imatrix& Derivative_lattice::get_smith_normal_form() const{
    
    return snf;
    
}

const imatrix& Derivative_lattice::get_row_matrix() const{
    
    return left;
    
}

const imatrix& Derivative_lattice::get_column_matrix() const{
    
    return right;
    
}

const imatrix& Derivative_lattice::get_product_table() const{
    
    return product_table;
    
}

const int& Derivative_lattice::get_snf_index() const{
    
    return snf_index;
    
}

const vector2_s& Derivative_lattice::get_label() const{
    
    return label;
    
}

void Derivative_lattice::set_label(const vector2_s& label_input){

    label = label_input;

}

