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

	    Header file for labeling.cpp

******************************************************************************/

#ifndef __LABELING
#define __LABELING

#include <vector>
#include <set>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "include/math.hpp"
#include "permutation.h"
#include "permutation_atom.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<short> > vector2_s;

class Labeling{

    vector2_s labeling;

    void combine_label(std::vector<short>& labeling_one_dim, 
        const std::vector<short>& labeling_single_lattice,
        const int& n_atom_labeling, const int& n_atom_labeling_single);

    void n_type_check(const std::vector<short>& label_candidate, 
        const std::vector<int>& n_type, const int& n_atom);

    int unimodular_check(const dmatrix& matrix, const double& symprec);

    vector2_i exchange_type(const std::vector<int>& n_type);

    public: 

    Labeling();
    ~Labeling();

    void all_labeling(const std::vector<int>& n_atom, 
            const std::vector<int>& n_type, 
            const std::vector<int>& n_atoms_derivative);

    void permutation_check(const vector2_s& label_candidate, 
            const imatrix& permutation_table);

    void unimodular_sym_check(
            const imatrix& hnf, const imatrix& right,
            const std::vector<imatrix>& rotate_matrix_array,
            const std::vector<dvector>& trans_vector_array,
            std::vector<imatrix>& rotate_matrix_candidate,
            std::vector<dvector>& trans_vector_candidate,
            const double& symprec);

    void exchange_check(
            const vector2_s& label_original, 
            const std::vector<int>& n_type,
            const imatrix& permutation);

    void superlattice_check(const vector2_s& label_candidate,
            const dmatrix& permutation_table_translation);

    const vector2_s& get_label() const;

};

#endif
