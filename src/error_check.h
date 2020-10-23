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

	    Header file for error_check.cpp

******************************************************************************/

#ifndef __ERROR_CHECK
#define __ERROR_CHECK

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

class Error_check{

	public: 

	Error_check();
	~Error_check();

    void check_derivative(const std::vector<int>& n_atoms_primitive, 
            const std::vector<int>& n_type);

    void check_correlation(const dmatrix& axis_primitive,
            const dmatrix& axis_change, const dmatrix& axis,
            const std::vector<int>& n_atoms_primitive,
            const std::vector<int>& n_atoms,
            const std::vector<int>& n_type,
            const std::vector<int>& spin_array,
            const double& symprec);

    void check_lsf_files(const dvector& energy, const dmatrix& correlation);

};

#endif
