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
        Fifth Floor, Boston, MA 02110-1301, USA, 
        or see http://www.gnu.org/copyleft/gpl.txt 

	    Program for searching for derivative structures

 **************************************************************************/

#include <iostream>
#include <vector>

#include "input.h"
#include "output.h"
#include "sym.h"

#include "hermite.h"
#include "smith.h"
#include "derivative_lattice.h"
#include "permutation_atom.h"
#include "labeling.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<short> > vector2_s; 
typedef std::vector<std::vector<std::vector<short> > > vector3_s; 

typedef std::vector<std::vector<int> > vector2_i; 
typedef std::vector<std::vector<std::vector<int> > > vector3_i; 

int main(){

    // reading input files
    Input ip;
    ip.input_derivative();

    const int n_lattice(ip.get_n_lattice());
    const std::vector<int> n_type(ip.get_n_type());
    const std::vector<int> n_atoms_derivative(ip.get_n_atoms_derivative());

    const dmatrix axis(ip.get_axis_primitive());
    const std::vector<int> n_atoms_primitive(ip.get_n_atoms_primitive());
    const std::vector<dvector> position(ip.get_position_primitive());
    const std::vector<int> type(ip.get_type_primitive());
    const double symprec(ip.get_symprec());

    // calculation of symmetric operations
    Sym symmetry(axis, position, type, symprec);
    const std::vector<imatrix> 
        rotate_matrix_array(symmetry.get_rotate_matrix());
    const std::vector<dvector> 
        trans_vector_array(symmetry.get_trans_vector());
    std::cout << "  number of symmetry operations = " 
        << rotate_matrix_array.size() << std::endl;

    // positions in n_sublattice (size of n_type)
    const std::vector<dvector> 
        position_simulation(ip.get_position_simulation());
    const int n_atom_primitive(position_simulation.size());
    //const int n_atom_label(n_atom_primitive * n_lattice);

    const std::vector<int> n_atoms(ip.get_n_atoms());

    // unique HNF
    Hermite herm(n_lattice, rotate_matrix_array, symprec);
    const std::vector<imatrix> hermite_array 
        (herm.get_hermite_normal_form_array());

    std::cout << "  number of unique superlattices (HNF) = "
        << hermite_array.size() << std::endl;

    // unique SNF and preparing Derivative_lattice
    std::vector<Derivative_lattice> derivative_lattice_array;
    std::vector<imatrix> smith_unique, permutation_unique;
    for (int i = 0; i < hermite_array.size(); ++i){
        Smith sm(hermite_array[i]);
        const imatrix snf(sm.get_smith_normal_form());
        const imatrix left(sm.get_row_matrix());
        const imatrix right(sm.get_column_matrix());
        Permutation_atom perm_atom;
        imatrix permutation 
            = perm_atom.permutation_lattice_translation_multilattice
            (snf, n_atom_primitive);

        int snf_index;
        if (i == 0){
            snf_index = 0;
            smith_unique.push_back(snf);
            permutation_unique.push_back(permutation);
        }
        else {
            for (int j = 0; j < smith_unique.size(); ++j){
                if (ublas::norm_inf(snf - smith_unique[j]) == 0){
                    snf_index = j;
                    break;
                }
                else if (j == smith_unique.size() - 1){
                    snf_index = smith_unique.size();
                    smith_unique.push_back(snf);
                    permutation_unique.push_back(permutation);
                }
            }
        }
        vector2_s label_empty;
        Derivative_lattice deriv(hermite_array[i], snf, left, right, 
                permutation, snf_index, label_empty);
        derivative_lattice_array.push_back(deriv);
    }

    std::cout << "  number of SNF = " << smith_unique.size() << std::endl;

    /*  Algorithm for finding unique labelings after obtaining HNF, 
        SNF and permutations by lattice translations.

    1. all labelings: in binary system, 2^n configurations (n sites)
        For a fixed composition, only configurations 
        which have the composition are selected.
        For a whole range of compositions, duplicates 
        by exchange operations are removed.
        This step is common for all HNF and SNF.

    2. Duplicates by lattice translations and superlattices which 
        can be expressed by using a smaller lattice are removed. 
        This step depends only on SNF.

    3. Duplicates by symmetric operations of primitive lattice are removed.
        Unimodular symmetric operations are only used for the duplicate check.
        The unimodular operations are dependent on HNF.
        Rotated labelings are examined by lattice translations.
    */
    

    // unique labeling for each unique SNF
    Labeling label;
    // all_labeling cannot work when n_atoms = 81 (Kuwabara)
    label.all_labeling(n_atoms, n_type, n_atoms_derivative);
    vector2_s labeling_init = label.get_label();
/*
    for (int i = 0; i < labeling_init.size(); ++i){
        for (int j = 0; j < labeling_init[i].size(); ++j){
            std::cout << labeling_init[i][j];
        }
        std::cout << std::endl;
    }
*/
    // duplication check by a lattice translation
    vector3_s all_label;
    for (int i = 0; i < smith_unique.size(); ++i){
        Labeling label;
        label.permutation_check(labeling_init, permutation_unique[i]);
        vector2_s labeling = label.get_label();
        all_label.push_back(labeling);
    }

    for (int i = 0; i < derivative_lattice_array.size(); ++i){
        int snf_index = derivative_lattice_array[i].get_snf_index();
        derivative_lattice_array[i].set_label(all_label[snf_index]);
    }
    all_label.clear();

    #ifdef _OPENMP
    #pragma omp parallel for 
    #endif
    for (int i = 0; i < derivative_lattice_array.size(); ++i){
        std::cout << "-----" << i << "--------" << std::endl;
        imatrix hnf = derivative_lattice_array[i].get_hermite_normal_form();
        imatrix snf = derivative_lattice_array[i].get_smith_normal_form();
        imatrix left = derivative_lattice_array[i].get_row_matrix();
        imatrix right = derivative_lattice_array[i].get_column_matrix();
        vector2_s label_original = derivative_lattice_array[i].get_label();
        std::cout << " HNF " << " = " << hnf << std::endl;
        std::cout << " SNF " << " = " << snf << std::endl;

       std::vector<imatrix> rotate_matrix_candidate;
        std::vector<dvector> trans_vector_candidate;

        Labeling label;
        label.unimodular_sym_check(hnf, right, rotate_matrix_array, 
            trans_vector_array, rotate_matrix_candidate, 
            trans_vector_candidate, symprec);

        Permutation_atom perm_atom;
        imatrix permutation_table = perm_atom.permutation_rotation
            (snf, left, rotate_matrix_candidate, trans_vector_candidate,
             position_simulation, symprec);


        label.permutation_check(label_original, permutation_table);
        vector2_s labeling = label.get_label();

        imatrix permutation_translation 
            = derivative_lattice_array[i].get_product_table();
        imatrix permutation = perm_atom.combine_table(permutation_table, 
            permutation_translation, snf);

        if (n_atoms_derivative.size() == 0) {
            label.exchange_check(labeling, n_type, permutation);
        }

        label.superlattice_check(label.get_label(), permutation_translation);

        labeling = label.get_label();
        derivative_lattice_array[i].set_label(label.get_label());
    }

    // output information of derivative structures
    Output op;
    op.output_derivative_out
        (derivative_lattice_array, axis, n_atoms_primitive, 
            position, n_type);

    return 0;

}
