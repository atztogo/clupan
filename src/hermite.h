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

	    Header file for hermite.cpp

******************************************************************************/

#ifndef __HERMITE
#define __HERMITE

#include <iostream>
#include <vector>
#include <set>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<int> ivector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

class Hermite{

    std::vector<imatrix> hermite_array;

    void candidate_hermite_normal_matrix(int n_lattice);
    void unique_hnf(const std::vector<imatrix>& rotate_matrix_array,
            const double& symprec);
    int unimodular_check(const dmatrix& matrix, const double& symprec);

    public: 

    Hermite(const int& n_lattice,
        const std::vector<imatrix>& rotate_matrix_array, const double& symprec);
    ~Hermite();

    const std::vector<imatrix>& get_hermite_normal_form_array() const;

};

#endif
