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

	    Class for obtaining convex hull of correlation functions
	
*****************************************************************************/

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "input.h"

namespace ublas = boost::numeric::ublas;
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
        std::cerr << " Specify name of correlation (option -f)." << std::endl;
        exit(8);
    }

    Input ip;
    ip.read_correlation(file_name);
    dmatrix correlation = ip.get_correlation();

    std::cout.precision(15);
    for (int i = 0; i < correlation.size1(); ++i){
        for (int j = 0; j < correlation.size2(); ++j){
//            double tmp = correlation(i,j) * pow(10, 15);
 ////           tmp = (double)(long)(tmp) * pow(10, -15) ;
//            if (fabs(tmp - 0.25) < 1e-15) tmp = 0.25;
//            else if (fabs(tmp - 0.5) < 1e-15) tmp = 0.5;
//            else if (fabs(tmp - 0.75) < 1e-15) tmp = 0.75;
//            else if (fabs(tmp - 1.0) < 1e-15) tmp = 1.0;
//            else if (fabs(tmp + 0.25) < 1e-15) tmp = -0.25;
//            else if (fabs(tmp + 0.5) < 1e-15) tmp = -0.5;
//            else if (fabs(tmp + 0.75) < 1e-15) tmp = -0.75;
//            else if (fabs(tmp + 1.0) < 1e-15) tmp = -1.0;
//            std::cout << tmp << " ";
            std::cout << round(correlation(i,j)*108216108000) << " ";
        }
        std::cout << std::endl;
    }


    return 0;

}
