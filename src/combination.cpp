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
        
        Class for combination

***************************************************************************/

#include "combination.h"

Combination::Combination(){}
Combination::~Combination(){}

imatrix Combination::get_combination(int n, int r){

    // size1 must be equal to comb(n, r)

    long int n_combination = n_comb(n, r);
    imatrix combination_output(n_combination, r);

    gsl_combination * c;
    c = gsl_combination_calloc (n, r);

    int n_count(0);
    do
    {
        for  (int i = 0; i < r; ++i){
            combination_output(n_count, i) = gsl_combination_get(c, i);
        }
        ++n_count;
    }
    while (gsl_combination_next (c) == GSL_SUCCESS);
    gsl_combination_free (c);

    return combination_output;

}

int Combination::n_comb(int n, int r){

    if ( n == r ) return 1;
    if ( r == 0 ) return 1;
    return n_comb(n-1, r-1) + n_comb(n-1, r);

}
