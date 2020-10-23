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

    Class for calculating electrostatic energy using Ewald method

 ******************************************************************************/

#include "ewald.h"

Ewald::Ewald(){}

Ewald::Ewald(const dmatrix& axis_input, 
        const std::vector<dvector>& position_input,
        const double& volume_input, const dmatrix& reciprocal_axis_input,
        const std::vector<double>& charge_input,
        const double& eta_input, const double& rmax_input, 
        const double& gmax_input):
    axis(axis_input), position(position_input), volume(volume_input),
    reciprocal_axis(reciprocal_axis_input),
    charge(charge_input), eta(eta_input), rmax(rmax_input), gmax(gmax_input)
{

    set_r_coordination();
    set_g_vector();

    for (int i = 0; i < position.size(); ++i){
        dvector pos = ublas::prod(axis, position[i]);
        position_cartesian.push_back(pos);
    }

    const double pi(3.14159265358979323846);

    g_energy_constant.resize(g_vector_array.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < g_vector_array.size(); ++i){
        dvector g_vector = g_vector_array[i];
        double square_g = ublas::inner_prod(g_vector, g_vector);
        double constant = 2 * pi / volume 
            * exp(- square_g / (4 * eta)) / square_g;
//        g_energy_constant.push_back(constant);
        g_energy_constant[i] = constant;
    }

}

Ewald::~Ewald(){}

void Ewald::calc_energy(){

    const std::complex<double> imag(0,1);
    const double pi(3.14159265358979323846);
    const double correction(energy_correction());

    double pot_r_sum(0.0), pot_g_sum(0.0), pot_s_sum(0.0);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:pot_r_sum, pot_g_sum, pot_s_sum)
#endif
    for (int i = 0; i < position_cartesian.size(); ++i){
        double pot_r(0.0);
        for (int j = 0; j < r_coordination[i].size(); ++j){
            int atom = r_coordination[i][j].first;
            double function = r_coordination[i][j].second;
            pot_r += charge[atom] * function;
        }

        std::complex<double> pot_g(0.0,0.0);
        for (int j = 0; j < position_cartesian.size(); ++j){
            dvector diff_pos = position_cartesian[j] - position_cartesian[i];
            for (int g = 0; g < g_vector_array.size(); ++g){
                double g_prod_r 
                    = ublas::inner_prod(g_vector_array[g], diff_pos);
                pot_g += charge[j] * g_energy_constant[g] 
                    * exp (imag * g_prod_r);
            }
        }

        double pot_r_tmp, pot_g_tmp, pot_s_tmp;
        pot_r_tmp = charge[i] * pot_r * 0.5;
        pot_g_tmp = charge[i] * pot_g.real();
        pot_s_tmp = - pow (charge[i], 2) * sqrt(eta/pi);
        pot_r_sum += pot_r_tmp;
        pot_g_sum += pot_g_tmp;
        pot_s_sum += pot_s_tmp;
    }

    pot_r_output = pot_r_sum * correction;
    pot_g_output = pot_g_sum * correction;
    pot_s_output = pot_s_sum * correction;
    pot_all = ( pot_r_sum + pot_g_sum + pot_s_sum ) * correction;
}

double Ewald::calc_mc_energy(const std::vector<double>& charge, 
    const std::vector<std::complex<double> >& charge_reciprocal,
    const int& atom1, const int& atom2){

    const double correction(energy_correction());
    const double diff_charge(charge[atom2] - charge[atom1]);
    const double diff_charge_square(pow (diff_charge, 2.0));

    double diff_pot_r(0.0);
    int r_size1 = r_coordination[atom1].size();
    int r_size2 = r_coordination[atom2].size();
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+:diff_pot_r)
#endif
    for (int i = 0; i < r_size1; ++i){
        int j = r_coordination[atom1][i].first;
        double function = r_coordination[atom1][i].second;
        if (j == atom1){
            diff_pot_r += diff_charge_square * function;
        }
        else if (j != atom2){
            diff_pot_r += diff_charge * charge[j] * function;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for reduction(+:diff_pot_r)
#endif
    for (int i = 0; i < r_size2; ++i){
        int j = r_coordination[atom2][i].first;
        double function = r_coordination[atom2][i].second;
        if (j == atom2){
            diff_pot_r -= diff_charge_square * function;
        }
        else if (j != atom1){
            diff_pot_r -= diff_charge * charge[j] * function;
        }
    }
    diff_pot_r *= correction;

//    std::complex<double> diff_pot_g(0.0,0.0);
    double diff_pot_g;
    const std::complex<double> imag(0,1);

#ifdef _OPENMP
#pragma omp parallel for reduction(+:diff_pot_g)
#endif
    for (int g = 0; g < g_vector_array.size(); ++g){
        dvector g_vector = g_vector_array[g];
        std::complex<double> pot1, pot2, pot3;
        double g_prod_r1 
            = ublas::inner_prod(g_vector, position_cartesian[atom1]);
        double g_prod_r2 
            = ublas::inner_prod(g_vector, position_cartesian[atom2]);
        double g_prod_r_diff = ublas::inner_prod
            (g_vector, position_cartesian[atom2]-position_cartesian[atom1]);
        pot1 = g_energy_constant[g] * charge_reciprocal[g]
            * (exp (-imag * g_prod_r1) - exp (-imag * g_prod_r2)) * diff_charge;
        pot2 = g_energy_constant[g] * pow(diff_charge, 2.0);
        pot3 = g_energy_constant[g] * pow(diff_charge, 2.0) 
            * exp (-imag * g_prod_r_diff);
        diff_pot_g += pot1.real() + pot2.real() - pot3.real();
    }
    diff_pot_g *=  2.0 * correction;

    //double diff_energy = diff_pot_r + diff_pot_g.real();
    double diff_energy = diff_pot_r + diff_pot_g;

    return diff_energy;
}

std::vector<std::complex<double> > Ewald::initial_calc_charge_reciprocal
(const std::vector<double>& charge){

    const std::complex<double> imag(0,1);

    std::vector<std::complex<double> > charge_reciprocal;
    charge_reciprocal.resize(g_vector_array.size());

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int g = 0; g < g_vector_array.size(); ++g){
        dvector g_vector(g_vector_array[g]);
        std::complex<double> charge_g(0.0, 0.0);
        for (int i = 0; i < position_cartesian.size(); ++i){
            double g_prod_r 
                = ublas::inner_prod(g_vector, position_cartesian[i]);
            charge_g += charge[i] * exp (imag * g_prod_r);
        }
        charge_reciprocal[g] = charge_g;
    }

    return charge_reciprocal;

}

void Ewald::calc_diff_charge_reciprocal(
        const std::vector<double>& charge,
        std::vector<std::complex<double> >& charge_reciprocal,
        const int& atom1, const int& atom2){

    const std::complex<double> imag(0,1);
    double diff_charge = charge[atom2] - charge[atom1];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int g = 0; g < charge_reciprocal.size() ; ++g){
        dvector g_vector(g_vector_array[g]);
        double g_prod_r1 
            = ublas::inner_prod(g_vector, position_cartesian[atom1]);
        double g_prod_r2 
            = ublas::inner_prod(g_vector, position_cartesian[atom2]);
        charge_reciprocal[g] += diff_charge
            * (exp (imag * g_prod_r1) - exp (imag * g_prod_r2));
    }
}


void Ewald::set_r_coordination(){

    dvector axis_length(3);
    for (int i = 0; i < 3; ++i){
        ublas::matrix_column<dmatrix> axis_component(axis, i);
        axis_length(i)
            = sqrt(ublas::inner_prod(axis_component, axis_component));
    }

    // revision is needed /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    dvector cell_expand(3);
    for (int i = 0; i < 3; ++i){
        double length(axis_length(i));
        cell_expand(i) = 0;
        while (rmax > length){
            length += axis_length(i);
            cell_expand(i) += 1;
        }
    }

    std::vector<dvector> lattice_tmp;
    for (int n1 = -3*cell_expand(0)-2; n1 < 3*cell_expand(0)+3; ++n1){
        for (int n2 = -3*cell_expand(1)-2; n2 < 3*cell_expand(1)+3; ++n2){
            for (int n3 = -3*cell_expand(2)-2; n3 < 3*cell_expand(2)+3; ++n3){
                dvector trans(3), trans_test(3);
                trans(0) = double(n1);
                trans(1) = double(n2);
                trans(2) = double(n3);
                dvector diff_cartesian = ublas::prod(axis, trans);
                double distance = sqrt(ublas::inner_prod
                        (diff_cartesian, diff_cartesian));
                if (distance < rmax){
                    lattice_tmp.push_back(trans);
                }
            }
        }
    }

    std::vector<dvector> lattice_translation;
    lattice_translation.push_back(lattice_tmp[0]);
    for (int i = 0; i < lattice_tmp.size(); ++i){
        for (int x = -1; x < 2; ++x){
            for (int y = -1; y < 2; ++y){
                for (int z = -1; z < 2; ++z){
                    dvector add(3), candidate;
                    add(0) = double(x);
                    add(1) = double(y);
                    add(2) = double(z);
                    candidate = lattice_tmp[i] + add;
                    for (int j = 0; j < lattice_translation.size(); ++j){
                        if (lattice_translation[j](0) == candidate(0)
                        and lattice_translation[j](1) == candidate(1)
                        and lattice_translation[j](2) == candidate(2)){
                            break;
                        }
                        else if (j==lattice_translation.size()-1){
                            lattice_translation.push_back(candidate);
                        }
                    }
                }
            }
        }
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    r_coordination.resize(position.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < position.size(); ++i){
        std::vector<std::pair<int, double> > atom_coordination;
        for (int j = 0; j < position.size(); ++j){
            dvector diff_pos = position[j] - position[i];
            for (int k = 0; k < lattice_translation.size(); ++k){
                dvector diff_cartesian 
                    = ublas::prod(axis, diff_pos + lattice_translation[k]);
                double distance = sqrt(ublas::inner_prod
                        (diff_cartesian, diff_cartesian));
                if (distance < rmax and fabs(distance) > 1e-10){
                    std::pair<int, double> tmp_pair;
                    tmp_pair.first = j;
                    tmp_pair.second = erfc(sqrt(eta) * distance) / distance;
                    atom_coordination.push_back(tmp_pair);
                }
            }
        }
//        r_coordination.push_back(atom_coordination);
        r_coordination[i] = atom_coordination;
    }
}

void Ewald::set_g_vector(){

    dvector reciprocal_axis_length(3);
    for (int i = 0; i < 3; ++i){
        ublas::matrix_column<dmatrix> axis_component(reciprocal_axis, i);
        reciprocal_axis_length(i)
            = sqrt(ublas::inner_prod(axis_component, axis_component));
    }

    std::vector<int> g_cell_expand;
    for (int i = 0; i < 3; ++i){
        double length(reciprocal_axis_length(i));
        int expand(1);
        while (gmax > length){
            length += reciprocal_axis_length(i);
            ++expand;
        }
        g_cell_expand.push_back(expand);
    }

    dvector g_vector(3), cartesian;
    double distance;
    for (int g1 =  -3*g_cell_expand[0]; g1 < 3*g_cell_expand[0]; ++g1){
        for (int g2 =  -3*g_cell_expand[1]; g2 < 3*g_cell_expand[1]; ++g2){
            for (int g3 =  -3*g_cell_expand[2]; g3 < 3*g_cell_expand[2]; ++g3){
                if (g1 != 0 or g2 != 0 or g3 != 0){
                    g_vector(0) = g1;
                    g_vector(1) = g2;
                    g_vector(2) = g3;
                    cartesian = ublas::prod(reciprocal_axis, g_vector);
                    distance = sqrt(ublas::inner_prod(cartesian, cartesian));
                    if (distance < gmax){
                        g_vector_array.push_back(cartesian);
                    }
                }
            }
        }
    }

}

double Ewald::energy_correction(){

    double energy_correction;
    const double ev(1.602177e-19);
    const double coulomb_electron(1.602176462e-19);
    const double coulomb_constant(8.9877742e9);
    energy_correction = coulomb_constant / ev 
        * 1e10 * pow(coulomb_electron, 2);

    return energy_correction;

}

const double& Ewald::get_real_energy() const{
    return pot_r_output;
}
const double& Ewald::get_reciprocal_energy() const{
    return pot_g_output;
}
const double& Ewald::get_self_energy() const{
    return pot_s_output;
}
const double& Ewald::get_energy() const{
    return pot_all;
}

