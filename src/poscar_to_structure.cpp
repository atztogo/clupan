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

	    Class for converting POSCAR to structure
	
*****************************************************************************/

#include <iostream>
#include <vector>

#include "getopt.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_poscar.h"
#include "output_structure.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

extern char *optarg;
extern int optind, opterr, optopt;

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
        std::cerr << " Specify name of POSCAR (option -f)." << std::endl;
        exit(8);
    }

    Parse_poscar parse(file_name);
    dmatrix axis = parse.get_axis();
    std::vector<dvector> position = parse.get_position();
    std::vector<int> n_atoms = parse.get_num_atoms();

    Output_structure out;
    out.output_from_structure(axis, n_atoms, position);

    return 0;

}
