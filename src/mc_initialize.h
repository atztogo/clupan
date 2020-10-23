/****************************************************************************

        Copyright (C) 2013 Atsuto Seko
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


	    Header file for mc_initialize.cpp
		
****************************************************************************/

#ifndef __MC_INITILIZE
#define __MC_INITILIZE

#include <vector>
#include <set>
#include <string>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "input.h"
#include "structure.h"
#include "sym.h"
#include "parse_structure.h"
#include "correlation.h"
#include "ewald.h"
#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<std::pair<ivector, int> > > vector2_pair;
typedef std::vector<std::vector<double> > vector2_d;
typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<std::vector<int> > > vector3_i;
typedef std::vector<std::vector<std::vector<std::vector<int> > > > vector4_i;

class MC_initialize{

    double pi;

    int n_step_ann, n_step_eqv;
   std::vector<int> spin;

    dvector eci, eci_without_empty;
    std::vector<dmatrix> cluster_function_array;
    std::vector<int> cluster_index;

    double energy;
    dvector correlation_array, n_cluster_array;
    std::vector<int> sublattice;
    dvector disorder_cf;

    vector3_i n_cluster_coordination;
    vector4_i coordination;

    // crystal structure
    dmatrix axis;
    std::vector<int> n_atoms; 
    std::vector<dvector> position;

    // Electrostatic energy calculation
    std::vector<double> charge;
    Ewald ewald_obj;

    void sort_eci_cluster
        (const Input& ip, vector2_pair& unique_cluster_array);

    std::vector<double> get_complete_comp_array
        (const std::vector<double>& comp_array, 
         const std::vector<int>& n_type);

    std::vector<int> init_structure_rand
        (const std::vector<int>& n_atoms, 
         const std::vector<int>& n_type,
         const std::vector<double>& complete_comp_array, 
         const std::vector<int>& spin_array,
         const double& symprec);

    std::vector<int> init_structure_file
        (const std::string& init_structure, 
         const std::vector<dvector>& position,
         const imatrix& supercell_matrix,
         const int& n_atom_simulation, 
         const std::vector<int>& spin_array,
         const double& symprec);

    void set_coordination_binary(std::vector<imatrix>& permutation_array);
    void set_coordination(std::vector<imatrix>& permutation_array);

    dvector calc_disorder_cf
        (const std::vector<double>& complete_comp_array,
         const vector2_pair& unique_cluster_array,
         const std::vector<int>& n_atoms_primitive,
         const std::vector<int>& n_type,
         const std::vector<int>& spin_array);

    std::vector<double> structure_to_charge
        (const std::vector<int>& spin, const std::vector<int>& spin_array, 
         const std::vector<int>& n_type, const std::vector<int>& n_atoms,
         const std::vector<double>& charge_array);

    public: 

    MC_initialize(const Input& ip);
    ~MC_initialize();

    // structure
    const dmatrix& get_axis() const;
    const std::vector<int>& get_n_atoms() const;
    const std::vector<dvector>& get_position() const;

    // size = number of atoms
    const std::vector<int>& get_spin() const;
    const std::vector<int>& get_sublattice() const;

    // input information
    const int& get_n_step_ann() const;
    const int& get_n_step_eqv() const;

    // for sqs
    const dvector& get_disorder_cf() const;

    // for energy and c.f. calculations
    const dvector& get_eci_without_empty() const;
    const std::vector<int>& get_cluster_index() const;
    const std::vector<dmatrix>& get_cluster_function_array() const;
    const dvector& get_n_cluster_array() const;
    const vector4_i& get_coordination() const;
    const vector3_i& get_n_cluster_coordination() const;

    // energy and c.f. of initial structure
    const double& get_energy() const;
    const dvector& get_correlation_array() const;

    // for calculation of electrostatic energy
    const std::vector<double>& get_charge() const;
    const Ewald& get_ewald_obj() const;
};

#endif
