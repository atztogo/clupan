/***************************************************************************

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
	
    	Header file for combination.cpp

***************************************************************************/

#ifndef __COMBINATION
#define __COMBINATION

#include <iostream>
#include <vector>

#include <gsl/gsl_combination.h>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<int> imatrix;

class Combination{

    int n_comb(int n, int r);

    public:

    Combination();
    ~Combination();

    imatrix get_combination(int n, int r);

};

#endif
