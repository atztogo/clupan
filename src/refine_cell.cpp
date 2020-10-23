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

	    Class for refining cell
	
*****************************************************************************/

#include <iostream>
#include <vector>

#include "getopt.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_structure.h"
#include "output_structure.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

extern char *optarg;
extern int optind, opterr, optopt;

extern "C" int spg_refine_cell
(double lattice[3][3], double position[][3], int types[], 
 const int num_atom, const double symprec );

int main(int argc, char *argv[]){

    const char *file_name = NULL;

    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"f:"))!=-1){
        switch( option ){
            case 'f' : file_name = optarg;
                       break;
        }
    }
    if (file_name == NULL){
        std::cerr << " Specify file name of structure (option -f)." 
            << std::endl;
        exit(8);
    }

    Parse_structure input(file_name);
    dmatrix axis(input.get_axis());
    std::vector<int> n_atoms(input.get_num_atoms());
    std::vector<dvector> position(input.get_position());
    std::vector<int> type = input.get_type();

    int natom = position.size();

    double symprec = 1e-5;

    double lattice[3][3];
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            lattice[i][j] = axis(i,j);
        }
    }

    int size = natom * 100000;
    double position2[size][3];
    int type2[size];

    for (int i = 0; i < natom; ++i){
        for (int j = 0; j < 3; ++j){
            position2[i][j] = position[i](j);
        }
        type2[i] = type[i];
    }

    int new_n_atom 
        = spg_refine_cell(lattice, position2, type2, natom, symprec);

    if (new_n_atom == 0){
        new_n_atom = natom;
    }

    dmatrix new_axis(3,3);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            new_axis(i,j) = lattice[i][j];
        }
    }

    std::vector<int> new_type;
    std::vector<dvector> new_position;
    for (int i = 0; i < new_n_atom; ++i){
        dvector pos(3);
        for (int j = 0; j < 3; ++j){
            pos(j) = position2[i][j];
        }
        new_position.push_back(pos);
        new_type.push_back(type2[i]);
    }

    Output_structure out;
    out.output_from_cell(new_axis, new_type, new_position);

    return 0;
}
