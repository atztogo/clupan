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

	    Main program for thermodynamic integration

*****************************************************************************/

#include <iostream>
#include <vector>
#include <iomanip>

#include "input.h"

typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<double> > vector2_d;

int main(int argc, char *argv[]){

    std::cout << " Definition of compsition and chemical potential "
        << "should be carefully checked." << std::endl;
    std::cout << " Definition of compsition and chemical potential "
        << "should be [mu(B) - mu(A)] comp(B)]" << std::endl;
    std::cout << " Source code should be modified using "
        << "Grand pot. = E-TS-[mu(B) - mu(A)* comp(B)" << std::endl;
    exit(8);

    const double k_B(0.86171 * pow (10,-4));

    Input ip;
    ip.input_ti(argc, argv);

    const double init_potential(ip.get_init_potential());
    const vector2_d mu_array(ip.get_ti_mu_array());
    const vector2_d temp_array(ip.get_ti_temp_array());
    const vector2_d comp_array(ip.get_ti_comp_array());
    const vector2_d energy_array(ip.get_ti_energy_array());

    for (int i = 0; i < mu_array.size() - 1; ++i){
        if ( mu_array[i][mu_array[i].size()-1] != mu_array[i+1][0]
        or temp_array[i][temp_array[i].size()-1] != temp_array[i+1][0]){
            std::cout << " End of path " << i + 1 
                << " is not consistent with the beginning of path " 
                << i + 2 << std::endl;
            exit(8);
        }
    }
    
    // converting temerature to beta = 1/(k_B * T)
    vector2_d beta_array;
    for (int i = 0; i < temp_array.size(); ++i){
        std::vector<double> tmp_vec;
        for (int j = 0; j < temp_array[i].size(); ++j){
            tmp_vec.push_back(1.0/k_B/temp_array[i][j]);
        }
        beta_array.push_back(tmp_vec);
    }

    // finding constant variable
    std::vector<int> switch_index;
    for (int i = 0; i < mu_array.size(); ++i){
        double sum_mu(0.0), sum_temp(0.0);
        for (int j = 0; j < mu_array[i].size() - 1; ++j){
            double diff_mu = mu_array[i][j+1] - mu_array[i][j];
            double diff_temp = temp_array[i][j+1] - temp_array[i][j];
            sum_mu += diff_mu;
            sum_temp += diff_temp;
        }
        if (fabs(sum_mu) > 1e-8 and fabs(sum_temp) > 1e-8){
            std::cout 
                << "Either mu or temperature must be constant in path "
                << i + 1 << " !" << std::endl;
            exit(8);
        }
        else if (fabs(sum_mu) > 1e-8 and fabs(sum_temp) < 1e-8){
            switch_index.push_back(0);
        }
        else {
            switch_index.push_back(1);
        }
    }

    // finding phase boundary
    int index(0);
    vector2_i phase_index;
    double boundary_criterion(0.2);
    for (int i = 0; i < comp_array.size(); ++i){
        std::vector<int> tmp_array;
        tmp_array.push_back(index);
        for (int j = 1; j < comp_array[i].size(); ++j){
            double diff_comp = fabs(comp_array[i][j] - comp_array[i][j-1]);
            if (diff_comp > boundary_criterion){
                ++index;
            }
            tmp_array.push_back(index);
        }
        phase_index.push_back(tmp_array);
    }

    // calculating grand_potential
    vector2_d grand_potential;
    double pot(beta_array[0][0] * init_potential);
    for (int i = 0; i < switch_index.size(); ++i){
        std::vector<double> pot_vec;
        if (switch_index[i] == 0){
            double beta = beta_array[i][0];
            pot_vec.push_back(pot / beta);
            double mu_0, mu_1, comp_0, comp_1;
            for (int j = 1; j < mu_array[i].size(); ++j){
                mu_0 = mu_array[i][j-1];
                mu_1 = mu_array[i][j];
                comp_0 = comp_array[i][j-1];
                comp_1 = comp_array[i][j];
                pot -= 0.5 * (mu_1 - mu_0) * beta * (comp_0 + comp_1);
                pot_vec.push_back(pot / beta);
            }
        }
        else if (switch_index[i] == 1){
            double mu = mu_array[i][0];
            pot_vec.push_back(pot / beta_array[i][0]);
            double beta_0, beta_1, e_mux_0, e_mux_1;
            for (int j = 1; j < mu_array[i].size(); ++j){
                beta_0 = beta_array[i][j-1];
                beta_1 = beta_array[i][j];
                e_mux_0 = energy_array[i][j-1] - mu * comp_array[i][j-1];
                e_mux_1 = energy_array[i][j] - mu * comp_array[i][j];
                pot += 0.5 * (beta_1 - beta_0) * (e_mux_0 + e_mux_1);
                pot_vec.push_back(pot / beta_array[i][j]);
            }
        }
        grand_potential.push_back(pot_vec);
    }

    // converting grand potential to Helmholtz free energy
    vector2_d helmholtz_free_energy;
    for (int i = 0; i < grand_potential.size(); ++i){
        std::vector<double> free_energy_vec;
        for (int j = 0; j < grand_potential[i].size(); ++j){
            double free_energy = grand_potential[i][j] 
                + (mu_array[i][j] * comp_array[i][j]);
            free_energy_vec.push_back(free_energy);
        }
        helmholtz_free_energy.push_back(free_energy_vec);
    }
            
    std::cout << std::setw(14) << std::right << "# chemical pot." 
        << " " << std::setw(14) << std::right << "temperature"  
        << " " << std::setw(15) << std::right << "grand pot."  
        << " " << std::setw(14) << std::right << "composition"
        << " " << std::setw(15) << std::right << "Helmholtz F" 
        << std::endl;
    std::cout << "#" << std::setw(13) << std::right << "  " 
        << " " << std::setw(14) << std::right << "(K)"  
        << " " << std::setw(15) << std::right << "(eV/atom)"  
        << " " << std::setw(14) << std::right << " "
        << " " << std::setw(15) << std::right << "(eV/atom)" 
        << std::endl;
    std::cout.precision(10);
    for (int i = 0; i < helmholtz_free_energy.size(); ++i){
        std::cout << "# path " << i + 1 << std::endl;
        for (int j = 0; j < helmholtz_free_energy[i].size(); ++j){
            if (phase_index[i][j] == 0){
                std::cout << std::setw(14) << std::right << mu_array[i][j] 
                    << " " << std::setw(14) << std::right << temp_array[i][j]  
                    << " " << std::setw(15) << std::right 
                    << grand_potential[i][j]  
                    << " " << std::setw(14) << std::right << comp_array[i][j]
                    << " " << std::setw(15) << std::right 
                    << helmholtz_free_energy[i][j] << std::endl;
            }
        }
    }

	return 0;

}
