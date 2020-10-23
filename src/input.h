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


	    Header file for input.cpp
		
****************************************************************************/

#ifndef __INPUT
#define __INPUT

#include <vector>
#include <string>
#include <numeric>
#include <iostream>
#include <sstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_input.h"
#include "parse_structure.h"
#include "structure.h"
#include "error_check.h"
#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<int> ivector;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<double> > vector2_d;
typedef std::vector<std::vector<std::pair<ivector, int> > > vector2_pair;

extern char *optarg;
extern int optind, opterr, optopt;

class Input{

    double pi;

    int n_lattice;
    double symprec;
    std::vector<int> n_type, n_atoms_derivative;

    std::vector<imatrix> hnf_array, snf_array, left_array;
    vector2_i label_array;

    int n_sublattice, max_cluster;
    std::vector<double> trunc_distance;

    dmatrix axis_change;
    vector2_pair unique_cluster_array;
    std::vector<int> spin_array;
    std::vector<int> spin_structure;
    char *file_name;

    dvector energy;
    dmatrix correlation, weight;
    std::vector<int> cluster_index;
    std::vector<std::string> energy_structure_index;

    int nloop, n_cluster;
    std::vector<double> temp_array;
    std::vector<int> base_array, cluster_candidate;
    int n_pop, max_iteration, n_elite;
    double p_mating, p_mutation;

    dvector eci;
    int composition_cluster;
    std::vector<double> composition_definition;
    int error_check;

    std::string init_structure;
    std::vector<dmatrix> cluster_function_array;
    std::vector<int> cluster_number_array;
    std::vector<double> comp_array, mu;
    int n_step_ann, n_step_eqv, i_ewald;
    std::string simulation_type;
    double epsilon;
    vector2_i n_change_array;

    std::vector<double> charge, charge_array;
    double eta, rmax, gmax;

    std::vector<int> group_index_array;
    std::vector<std::string> structure_index_array;

    double init_potential;
    vector2_d ti_mu_array, ti_temp_array, ti_comp_array, ti_energy_array;

    std::string criterion;

    // crystal structure
    dmatrix axis, axis_primitive, reciprocal_axis;
    std::vector<int> n_atoms, n_atoms_primitive, type, type_primitive;
    std::vector<dvector> position, position_primitive, position_simulation;
    double volume;

    void read_cluster_out();
    void read_weight(const char* file_name);
    void read_eci(const char* file_name);
    void read_cluster_function_out();

    void set_position_simulation
        (const int& n_sublattice,
         const std::vector<dvector>& position_all, 
         const std::vector<int>& n_atoms);

    void get_temp_array(Parse_input& input, const int& i_temp);
    void read_ewald_input
        (Parse_input& input, const int& n_atom, const double& volume);

    void error_spin_array_size
        (const std::vector<int>& spin_array, const std::vector<int>& n_type);
    void error_comp_array_size
        (const std::vector<double>& comp_array,
         const std::vector<int>& spin_array, const std::vector<int>& n_type);
    std::string read_file_foption(int argc, char *argv[]);
    void set_weight(Parse_input& parse);

    public: 

    Input();
    ~Input();

    void input_derivative();
    void input_derivative_to_structure();
    void input_derivative_to_structure_fix_comp();
    void input_cluster();
    void input_correlation(int argc, char *argv[]);
    void input_lsf();
    void input_cv_casp(int argc, char *argv[]);
    void input_ga();
    void input_gss();
    void input_mc();
    void input_ewald(int argc, char *argv[]);
    void input_ti(int argc, char *argv[]);
    void input_sqs();

    void read_energy(const char* file_name);
    void read_correlation(const char* file_name);
    void read_group(const char* file_name);

    // for derivative
    const std::vector<int>& get_n_type() const;
    const int& get_n_lattice() const;
    const double& get_symprec() const;
    const std::vector<int>& get_n_atoms_derivative() const;

    // for derivative_to_structure
    const std::vector<imatrix>& get_hnf_array() const;
    const std::vector<imatrix>& get_snf_array() const;
    const std::vector<imatrix>& get_left_array() const;
    const vector2_i& get_label_array() const;

    // for cluster_search
    const int& get_max_cluster() const;
    const int& get_n_sublattice() const;
    const std::vector<double>& get_trunc_distance() const;

    // for correlation
    const dmatrix& get_axis_change() const;
    const vector2_pair& get_unique_cluster_array() const;
    const std::vector<int>& get_spin_array() const;
    const std::vector<int>& get_spin_structure() const;
    const char* get_file_name() const;

    // for lsf
    const dvector& get_energy() const;
    const dmatrix& get_correlation() const;
    const dmatrix& get_weight() const;
    const std::vector<std::string>& get_energy_structure_index() const;
    const std::vector<int>& get_cluster_index() const;

    // for ga
    const int& get_nloop() const;
    const int& get_n_cluster() const;
    const std::vector<double>& get_temp_array() const;
    const std::vector<int>& get_base_array() const;
    const std::vector<int>& get_cluster_candidate() const;
    const int& get_n_pop() const;
    const int& get_max_iteration() const;
    const int& get_n_elite() const;
    const double& get_p_mating() const;
    const double& get_p_mutation() const;

    const dvector& get_eci() const;
    const int& get_composition_cluster() const;
    const std::vector<double>& get_composition_definition() const;
    const int& get_error_check() const;
    const std::vector<dmatrix>& get_cluster_function_array() const;
    const std::vector<int>& get_cluster_number_array() const;
    const std::string& get_init_structure() const;
    const std::vector<double>& get_comp_array() const;
    const int& get_n_step_ann() const;
    const int& get_n_step_eqv() const;
    const std::string& get_simulation_type() const;
    const int& get_i_ewald() const;
    const double& get_epsilon() const;
    const std::vector<double>& get_mu() const;
    const std::vector<double>& get_charge() const;
    const std::vector<double>& get_charge_array() const;
    const vector2_i& get_n_change_array() const;
    const double& get_weight_factor() const;
    const double& get_accuracy() const;
    const double& get_eta() const;
    const double& get_rmax() const;
    const double& get_gmax() const;
    const std::vector<std::string>& get_structure_index_array() const;
    const std::vector<int>& get_group_index_array() const;
    const double& get_init_potential() const;
    const vector2_d& get_ti_mu_array() const;
    const vector2_d& get_ti_temp_array() const;
    const vector2_d& get_ti_comp_array() const;
    const vector2_d& get_ti_energy_array() const;
    const std::string& get_criterion() const;

    const dmatrix& get_axis() const;
    const std::vector<int>& get_n_atoms() const;
    const std::vector<dvector>& get_position() const;
    const std::vector<int>& get_type() const;
    const dmatrix& get_axis_primitive() const;
    const std::vector<int>& get_n_atoms_primitive() const;
    const std::vector<dvector>& get_position_primitive() const;
    const std::vector<int>& get_type_primitive() const;
    const std::vector<dvector>& get_position_simulation() const;
    const double& get_volume() const;
    const dmatrix& get_reciprocal_axis() const;


};

#endif
