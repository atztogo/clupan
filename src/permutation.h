/*****************************************************************************

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

	    Header file for permutation.cpp

******************************************************************************/

#ifndef __PERMUTATION
#define __PERMUTATION

#include <iostream>
#include <vector>
#include <cmath>

#include <gsl/gsl_permutation.h>

typedef std::vector<std::vector<int> > vector2_i;

class Permutation{

    vector2_i permutation;
    std::vector<short> permutation_one_dimension;

    int composition_check(short comb[ ], 
        const std::vector<int>& n_atoms_derivative);

    public: 

    Permutation();
    ~Permutation();

    void permutation_label(const int& n_atom, const int& n_type);
    void permutation_label_fix_composition
        (const std::vector<int>& n_atoms_derivative);

    void permutation_without_repetition(const int& n_type);
    void permutation_with_repetition(const int& n_atom, const int& n_type);
    void permutation_with_repetition_multilattice(
            const std::vector<int>& n_type_array);

    const vector2_i& get_permutation() const;
    const std::vector<short>& get_permutation_label() const;

};

#endif
