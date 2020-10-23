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

    	Header file for least_square.cpp

******************************************************************************/

#ifndef __LEAST_SQUARES
#define __LEAST_SQUARES

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>

//#include "include/math.hpp"

typedef boost::numeric::ublas::vector<double> dvector;
typedef boost::numeric::ublas::matrix<double> dmatrix;

namespace ublas = boost::numeric::ublas;

class Least_squares{

    dvector eci;
    double square_error;

    public: 

    Least_squares();
    ~Least_squares();

    int fit(const dvector& energy, const dmatrix& correlation, 
        const dmatrix& weight);
    const dvector& get_eci() const;
    const double& get_square_error();
};

#endif
