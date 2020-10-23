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

	
	    Class for checking errors

******************************************************************************/

#include "error_check.h"

Error_check::Error_check(){}

Error_check::~Error_check(){}

void Error_check::check_derivative(
        const std::vector<int>& n_atoms_primitive, 
        const std::vector<int>& n_type){

    if (n_type.size() > n_atoms_primitive.size()){
        std::cerr << " error : size of n_type is larger than size of n-atoms."
            << std::endl;
        exit(8);
    }
    
}

void Error_check::check_correlation(const dmatrix& axis_primitive, 
        const dmatrix& axis_change, const dmatrix& axis, 
        const std::vector<int>& n_atoms_primitive,
        const std::vector<int>& n_atoms, 
        const std::vector<int>& n_type, 
        const std::vector<int>& spin_array, 
        const double& symprec){

    // axis check
    dmatrix axis_primitive_prod = ublas::prod(axis_primitive, axis_change);
    if (ublas::norm_inf(axis - axis_primitive_prod) > symprec){
        std::cerr << " error : axis != axis_primitive * axis_change "
            << std::endl;
        exit(8);
    }
    // n_atoms check
    int n_all_atoms(0), n_all_atoms_primitive(0);
    for (int i = 0; i < n_atoms.size(); ++i){
        n_all_atoms += n_atoms[i];
    }
    for (int i = 0; i < n_atoms_primitive.size(); ++i){
        n_all_atoms_primitive += n_atoms_primitive[i];
    }
    int remainder = n_all_atoms % n_all_atoms_primitive;
    int expansion = n_all_atoms / n_all_atoms_primitive;
    double det = math::determinant(axis_change);
    if (remainder != 0){ 
        std::cerr << " n-atoms in structure cannot be ";
        std::cerr << "divided by n-atom in primitive cell." << std::endl;
        exit(8);
    }
    if (fabs(expansion - fabs(det)) > symprec){ 
        std::cerr << " cell sizes of structure and ";
        std::cerr << "primitive cell are different. " << std::endl;
        exit(8);
    }

    // spin array size check
    int sum(0);
    for (int i = 0; i < n_type.size(); ++i){
        sum += n_type[i];
    }
    if (sum != spin_array.size()){
        std::cerr << " sum of n_type in structure is different ";
        std::cerr << "from size of spin in correlation.in" << std::endl;
        exit(8);
    }    

}

void Error_check::check_lsf_files(const dvector& energy, 
    const dmatrix& correlation){

    if (correlation.size1() != energy.size()){
        std::cerr << "  number of structures in correlation = "
            << correlation.size1() << std::endl;
        std::cerr << "  number of structures in dft_energy = "
            << energy.size() << std::endl;
        exit(8);
    }
}

