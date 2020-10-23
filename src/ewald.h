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

	    Header file for ewald.cpp

******************************************************************************/

#ifndef __EWALD
#define __EWALD

#include <vector>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

class Ewald{

    dmatrix axis, reciprocal_axis;
    std::vector<dvector> position, position_cartesian;

    double volume;
    std::vector<double> charge;
    double eta, rmax, gmax;

    std::vector<std::vector<std::pair<int, double> > > r_coordination;
    std::vector<dvector> g_vector_array;
    std::vector<double> g_energy_constant;
    double pot_all, pot_r_output, pot_g_output, pot_s_output;

    void set_r_coordination();
    void set_g_vector();
    double energy_correction();

    public: 

    Ewald();
    Ewald(const dmatrix& axis_input, 
            const std::vector<dvector>& position_input, 
            const double& volume_input, 
            const dmatrix& reciprocal_axis_input,
            const std::vector<double>& charge_input,
            const double& eta_input, 
            const double& rmax_input, 
            const double& gmax_input);
    ~Ewald();

    void calc_energy();

    //void initialize_for_mc();
    double calc_mc_energy(const std::vector<double>& charge, 
            const std::vector<std::complex<double> >& charge_reciprocal,
            const int& atom1, const int& atom2);

    std::vector<std::complex<double> > initial_calc_charge_reciprocal
        (const std::vector<double>& charge);
    void calc_diff_charge_reciprocal(
        const std::vector<double>& charge,
        std::vector<std::complex<double> >& charge_reciprocal,
        const int& atom1, const int& atom2);

    const double& get_real_energy() const;
    const double& get_reciprocal_energy() const;
    const double& get_self_energy() const;
    const double& get_energy() const;

};

#endif
