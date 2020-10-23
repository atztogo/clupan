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


	    Header file for output_structure.cpp
		
****************************************************************************/

#ifndef __OUTPUT_STRUCTURE
#define __OUTPUT_STRUCTURE

#include <vector>
#include <iostream>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

class Output_structure{

    public: 

    Output_structure();
    ~Output_structure();

    void output_from_structure(const dmatrix& axis, 
        const std::vector<int>& n_atoms, 
        const std::vector<dvector>& position);

    void output_from_structure_to_file(std::ofstream& structure_out,
            const dmatrix& axis, const std::vector<int>& n_atoms, 
            const std::vector<dvector>& position);

    void output_from_cell(const dmatrix& axis, 
        const std::vector<int>& type, const std::vector<dvector>& position);

    void output_from_cell_to_file(std::ofstream& structure_out,
            const dmatrix& axis, const std::vector<int>& type, 
            const std::vector<int>& order, 
            const std::vector<dvector>& position);

};

#endif
