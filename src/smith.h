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

	    Header file for smith.cpp

******************************************************************************/

#ifndef __SMITH
#define __SMITH

#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<int> imatrix;

class Smith{

    imatrix snf, row_matrix, column_matrix, product_table;

    void calc_smith_normal_form(int size);
    void calc_product_table();
    imatrix permutation_matrix
        (const int& index1, const int& index2, const int& size);
    imatrix multiply_matrix
        (const int& index1, const int& multiplier, const int& size);
    imatrix addition_matrix
        (const int& index1, const int& index2, 
        const int& multiplier, const int& size);

    public: 

    Smith(const imatrix& matrix);
    ~Smith();

    const imatrix& get_smith_normal_form() const;
    const imatrix& get_row_matrix() const;
    const imatrix& get_column_matrix() const;
    const imatrix& get_product_table() const;

};

#endif
