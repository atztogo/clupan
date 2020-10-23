/************************************************************************

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
	

	    Header file for converting a label expression to another expression
		
************************************************************************/

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<int> imatrix;
typedef ublas::vector<int> ivector;

#ifndef __LABEL
#define __LABEL

template<typename T>
void sequential_number_to_lattice_expression(const T& sequential_number, 
    const imatrix& snf, ivector& lattice, int& atom_index){

    lattice.resize(3);

    int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);

    atom_index = int(sequential_number / n_lattice);
    int remainder = sequential_number % n_lattice;

    lattice(0) = remainder / (snf(1,1) * snf(2,2));
    lattice(1) = (remainder / snf(2,2) ) % snf(1,1);
    lattice(2) = remainder % snf(2,2);

}

template<typename T> 
void lattice_expression_to_sequential_number(T& lattice, 
    const int& atom_index, const imatrix& snf, int& sequential_number){

    int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);

    sequential_number = lattice(0) * snf(1,1) * snf(2,2) 
        + lattice(1) * snf(2,2) + lattice(2);

    sequential_number += n_lattice * atom_index;

}

#endif
