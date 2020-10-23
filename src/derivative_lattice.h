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

	    Header file for derivative_lattice.cpp

******************************************************************************/

#ifndef __DERIVATIVE_LATTICE
#define __DERIVATIVE_LATTICE

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;
typedef std::vector<std::vector<short> > vector2_s;

class Derivative_lattice{

    imatrix hnf, snf, left, right, product_table;
    int snf_index;
    vector2_s label;

    public: 

    Derivative_lattice(const imatrix& hnf_input, const imatrix& snf_input,
        const imatrix& left_input, const imatrix& right_input,
        const imatrix& product_table_input, const int& snf_index_input, 
        const vector2_s& label_input);
    ~Derivative_lattice();

    const imatrix& get_hermite_normal_form() const;
    const imatrix& get_smith_normal_form() const;
    const imatrix& get_row_matrix() const;
    const imatrix& get_column_matrix() const;
    const imatrix& get_product_table() const;
    const int& get_snf_index() const;
    const vector2_s& get_label() const;
    void set_label(const vector2_s& label_input);

};

#endif
