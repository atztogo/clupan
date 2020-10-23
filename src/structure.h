/****************************************************************************

        Copyright (C) 2013 Atsuto Seko
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


	    Header file for structure.cpp
		
****************************************************************************/

#ifndef __STRUCTURE
#define __STRUCTURE

#include <vector>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"
#include "include/label.hpp"
#include "smith.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::vector<int> ivector;
typedef ublas::matrix<double> dmatrix;
typedef ublas::matrix<int> imatrix;

class Structure{

    dmatrix axis;
    std::vector<dvector> position;
    std::vector<int> n_atoms;
    double volume;

    public: 
    Structure
        (const dmatrix& axis_input, 
         const std::vector<dvector>& position_input, 
         const std::vector<int>& n_atoms_input);
    ~Structure();

    const dmatrix& get_axis() const;
    const std::vector<int>& get_n_atoms() const;
    const std::vector<dvector>& get_position() const;
    double get_volume();
    dmatrix get_reciprocal_axis();
    int supercell(const imatrix& supercell_matrix);

};

#endif
