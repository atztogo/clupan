/****************************************************************************

        Copyright (C) 2007 Atsuto Seko
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


	    Header file for output.cpp
		
****************************************************************************/

#ifndef __OUTPUT
#define __OUTPUT

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "derivative_lattice.h"
#include "output_structure.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<int> ivector;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<short> > vector2_s;
typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<double> > vector2_d;
typedef std::vector<std::vector<dmatrix> > vector2_dmatrix;
typedef std::vector<std::vector<std::pair<ivector, int> > > vector2_pair;

class Output{

    void mc_spin_to_structure
        (const std::vector<int>& n_type, 
         const std::vector<int>& n_atoms, 
         std::vector<int>& spin, 
         std::vector<int>& spin_array);

    public: 

    Output();
    ~Output();

    void output_derivative_out(
            const std::vector<Derivative_lattice>& derivative_lattice_array,
            const dmatrix& axis, const std::vector<int>& num_atoms,
            const std::vector<dvector>& position, 
            const std::vector<int>& n_type);

    void output_derivative_structure(
            const std::vector<dmatrix>& new_axis_array, 
            const vector2_i& new_type_array,
            const std::vector<std::vector<dvector> >& new_position_array,
            const std::vector<dmatrix>& axis_change_array,
            const std::vector<int>& n_type);

    void output_cluster_out(
            const vector2_pair& all_unique_cluster_output, 
            const vector2_i& all_unique_cluster, 
            const dmatrix& axis_primitive, const std::vector<int>& n_atoms, 
            const std::vector<dvector>& position_primitive, 
            const int& n_sublattice, const dmatrix& distance_array);

    void output_correlation(const dvector& correlation_array, 
            const std::vector<int>& n_clusters,
            const std::vector<dmatrix>& cluster_function_array, 
            const std::vector<int>& cluster_number_array, 
            const char* structure_name,
            const std::vector<int>& n_type);

    void output_lsf(const std::vector<int>& cluster_index,
            const dvector& eci, const double& square_error, 
            const double& cv_score, const dmatrix& inverse, 
            const int& inverse_check);

    void output_cv_casp(const std::vector<double>& cv_casp, 
        const std::vector<int> n_structure_in_group, const double& cv_score);

    void output_ga(const vector2_i& pop_array, 
            const std::vector<double>& score_array, const int& iter);

    void output_gss(const std::vector<double>& composition_array, 
            const std::vector<double>& energy_array,
            const std::vector<double>& error_array,
            const double& rms_error,
            const std::vector<double>& gs_comp,
            const std::vector<double>& gs_energy,
            const std::vector<int>& gs_line);

    void output_mc
        (std::vector<int>& spin, const double& energy_fin, 
         const double& energy_ave, const double& specific_heat, 
         const dvector& correlation_ave, 
         const std::vector<double>& temp_array, 
         const dmatrix& axis, 
         const std::vector<int>& n_atoms, 
         const std::vector<dvector>& position, 
         const std::vector<int>& cluster_index, 
         const std::vector<int>& n_type, 
         std::vector<int>& spin_array);

    void output_sqs(const dvector& correlation_array, 
            const std::vector<int>& cluster_index, 
            const dmatrix& axis, const std::vector<int>& n_atoms, 
            const std::vector<dvector>& position, 
            const std::vector<int>& n_type, 
            std::vector<int>& spin, std::vector<int>& spin_array);

};

#endif
