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

    Main program for determining cell size for  Monte Carlo simulation

******************************************************************************/

#include <iostream>
#include <vector>

#include "parse_input.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<double> dmatrix;

int gcd(int m, int n);
int lcm(int m, int n);

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc, char *argv[]){

    const char* file_name_structure = NULL;
    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"f:"))!=-1){
        switch( option ){
            case 'f' : file_name_structure = optarg;
                       break;
        }
    }

    if (file_name_structure == NULL){
        std::cerr << " Specify file name of initial structure "
            << "for mc simulation (option -f)." << std::endl;
        exit(8);
    }

    dmatrix axis_change;
    Parse_input input(file_name_structure);
    input.assignNeed("axis_change", axis_change);

    int tmp = axis_change(2,2);
    for (int i = 0; i < 2; ++i){
        if (axis_change(2,i) != 0){
            tmp = lcm(axis_change(2,i), tmp);
        }
    }
    for (int i = 0; i < 3; ++i){
        if (axis_change(2,i) != 0){
            int test = tmp / axis_change(2,i);
            for (int j = 0; j < 3; ++j){
                axis_change(j,i) *= test;
            }
        }
    }

    tmp = axis_change(1,1);
    if (axis_change(1,0) != 0){
        tmp = lcm(axis_change(1,0), axis_change(1,1));
    }
    for (int i = 0; i < 3; ++i){
        if (axis_change(1,i) != 0){
            int test = tmp / axis_change(1,i);
            for (int j = 0; j < 3; ++j){
                axis_change(j,i) *= test;
            }
        }
    }

    std::cout << " " << file_name_structure << " can be expressed by ";
    for (int i = 0; i < 3; ++i){
        std::cout << axis_change(i,i) << " ";
    }
    std::cout << "expansion of primitive lattice. " << std::endl;

    return 0;

}

int gcd(int m, int n){
    int r;
    if (m < n) {
        int tmp = m;
        m = n;
        n = tmp;
    }
    if (n == 0) {
        return m;
    }
    if ((r = m % n)) {
        return gcd(n, r);
    }
    return n;
}

int lcm(int m, int n){
    return (m * n) / gcd(m, n);
}

