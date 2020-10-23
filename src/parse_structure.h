/****************************************************************************

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


	    Header file for parse_structure.cpp
		
****************************************************************************/

#ifndef __PARSE_STRUCTURE
#define __PARSE_STRUCTURE

#include <vector>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_input.h"
#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

class Parse_structure{

    dmatrix axis;
    std::vector<dvector> position;
    std::vector<int> n_atoms;
    std::vector<int> type;
    double volume;

    double modify_position(double x) const;

    public: 
    Parse_structure(const char *filename);
    ~Parse_structure();

    const dmatrix& get_axis() const;
    const std::vector<int>& get_num_atoms() const;
    const std::vector<dvector>& get_position() const;
    const std::vector<int>& get_type() const;
    double get_volume();
    dmatrix get_reciprocal_axis();

};

#endif
