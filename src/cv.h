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


	Header file for cv.cpp

*****************************************************************************/

#ifndef __CV
#define __CV

#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;

typedef boost::numeric::ublas::matrix<double> dmatrix;
typedef boost::numeric::ublas::vector<double> dvector;
typedef boost::numeric::ublas::vector<long double> ldvector;

class CV{

    double cv_score;
    std::vector<double> cv_casp;
    std::vector<int> n_structure_in_group;

    public: 

    CV();
    ~CV();

    double calc_cv_score(const dvector& energy, 
            const dmatrix& correlation, const dmatrix& weight, 
            const dvector& eci);

    std::vector<double> calc_cv_casp(const dvector& energy, 
            const dmatrix& correlation, const dvector& eci, 
            const std::vector<int>& group);
    const std::vector<int>& get_n_structure_in_group() const;

};

#endif
