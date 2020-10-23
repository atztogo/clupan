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

	    Main program for ground state search

*****************************************************************************/

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "input.h"
#include "output.h"

namespace ublas = boost::numeric::ublas;
typedef boost::numeric::ublas::vector<double> dvector;
typedef boost::numeric::ublas::matrix<double> dmatrix;

int main(){

    Input ip;
    ip.input_gss();
    const std::vector<int> cluster_index(ip.get_cluster_index());
    const dvector eci(ip.get_eci());
    const dmatrix correlation(ip.get_correlation());
    const int composition_cluster(ip.get_composition_cluster());
    const std::vector<double> composition_definition 
        (ip.get_composition_definition());

    // calculating CE energies of all structures
    std::vector<double> composition_array, energy_array;
    double composition, energy;
    for (int i = 0; i < correlation.size1(); ++i){
        composition = 0.0;
        for (int j = 0; j < composition_definition.size(); ++j){
            composition += composition_definition[j] 
                * pow(correlation(i, composition_cluster), j);
        }
        energy = 0.0;
        for (int j = 0; j < cluster_index.size(); ++j){
            energy += eci(j) * correlation(i, cluster_index[j]);
        }
        composition_array.push_back(composition);
        energy_array.push_back(energy);
    }

    // finding ground state structures
    std::set<double> comp_set;
    for (int i = 0; i < composition_array.size(); ++i){
        comp_set.insert(composition_array[i]);
    }

    std::vector<std::pair<double,int> > lowest_energy;
    std::set<double>::iterator it = comp_set.begin();
    while (it != comp_set.end()){
        std::map<double,int> all_energy_comp;
        for (int i = 0; i < composition_array.size(); ++i){
            if (*it == composition_array[i]){
                std::pair<double,int> tmp;
                tmp.first = energy_array[i];
                tmp.second = i;
                all_energy_comp.insert(tmp);
            }
        }
        lowest_energy.push_back(*all_energy_comp.begin());
        ++it; 
    }

    int index1, index2, index3;
    double energy1, energy2, energy3, comp1, comp2, comp3;
    double ratio2, ratio3;
    double energy_separation;
    std::vector<double> gs_comp, gs_energy;
    std::vector<int> gs_line;
    gs_comp.push_back(composition_array[lowest_energy[0].second]);
    gs_energy.push_back(lowest_energy[0].first);
    gs_line.push_back(lowest_energy[0].second + 1);
    for (int i = 1; i < lowest_energy.size() - 1; ++i){
        energy1 = lowest_energy[i].first;
        index1 = lowest_energy[i].second;
        comp1 = composition_array[index1];
        for (int j = 0; j < i; ++j){
            energy2 = lowest_energy[j].first;
            index2 = lowest_energy[j].second;
            comp2 = composition_array[index2];
            for (int k = i + 1; k < lowest_energy.size(); ++k){
                energy3 = lowest_energy[k].first;
                index3 = lowest_energy[k].second;
                comp3 = composition_array[index3];
                ratio2 = (comp3-comp1)/(comp3-comp2);
                ratio3 = (comp1-comp2)/(comp3-comp2);
                energy_separation = energy2 * ratio2 + energy3 * ratio3;
                if (energy1 > energy_separation - 1e-10){
                    goto loop;
                }
                else if (j == i - 1 and k == lowest_energy.size() - 1){
                    gs_comp.push_back(composition_array[index1]);
                    gs_energy.push_back(energy1);
                    gs_line.push_back(index1+1);
                }
            }
        }
        loop:;
    }
    gs_comp.push_back(composition_array[(*lowest_energy.rbegin()).second]);
    gs_energy.push_back((*lowest_energy.rbegin()).first);
    gs_line.push_back((*lowest_energy.rbegin()).second + 1);

    // calculating differences between CE and DFT energies
    std::vector<double> error_array;
    double rms_error(0.0);

    const int error_check(ip.get_error_check());
    if (error_check == 1){
        const dvector dft_energy(ip.get_energy());
        for (int i = 0; i < dft_energy.size(); ++i){
            double diff = energy_array[i] - dft_energy(i);
            error_array.push_back(diff);
            rms_error += pow (diff, 2.0);
        }
        rms_error /= double(dft_energy.size());
        rms_error = sqrt(rms_error);
    }

    Output op;
    op.output_gss(composition_array, energy_array, error_array, rms_error,
        gs_comp, gs_energy, gs_line);

	return 0;

}

