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

	    Class for generating output files
	
*****************************************************************************/

#include "output.h"

Output::Output(){}

Output::~Output(){}

void Output::output_derivative_out(
    const std::vector<Derivative_lattice>& derivative_lattice_array,
    const dmatrix& axis, const std::vector<int>& num_atoms, 
    const std::vector<dvector>& position, 
    const std::vector<int>& n_type){

    std::ofstream output;
    output.open("derivative.out",std::ios::out);
    output << " n_type = ";
    for (int i = 0; i < n_type.size(); ++i){
        output << n_type[i] << " ";
    }
    output << std::endl;
   // output << " n_sublattice = " << n_sublattice << std::endl;

    Output_structure os;
    os.output_from_structure_to_file(output, axis, num_atoms, position);

    int n_label(1);
    output << "# index labeling   HNF   SNF   left" << std::endl;
    for (int i = 0; i < derivative_lattice_array.size(); ++i){
        imatrix hnf = derivative_lattice_array[i].get_hermite_normal_form();
        imatrix snf = derivative_lattice_array[i].get_smith_normal_form();
        imatrix left = derivative_lattice_array[i].get_row_matrix();
        vector2_s label = derivative_lattice_array[i].get_label();

        for (int j = 0; j < label.size(); ++j){
            output << n_label << " ";
            for (int k = 0; k < label[j].size(); ++k){
                output << label[j][k];
            }
            output << " "  << hnf << " " << snf << " " << left << std::endl;
            ++n_label;
        }
    }
    output.close();

}

void Output::output_derivative_structure(
    const std::vector<dmatrix>& new_axis_array, 
    const vector2_i& new_type_array,
    const std::vector<std::vector<dvector> >& new_position_array,
    const std::vector<dmatrix>& axis_change_array, 
    const std::vector<int>& n_type){

    std::set<int> tmp;
    for (int i = 0; i < new_type_array[0].size(); ++i){
        tmp.insert(new_type_array[0][i]);
    }
    std::vector<int> order;
    order.resize(tmp.size());
    std::copy(tmp.begin(), tmp.end(), order.begin());

    for (int i = 0; i < new_axis_array.size(); ++i){
        dmatrix axis = new_axis_array[i];
        std::vector<int> type = new_type_array[i];
        std::vector<dvector> position = new_position_array[i];
        dmatrix axis_change = axis_change_array[i];

        std::string str("structure"), str2;
        std::stringstream ss;
        ss << i + 1;
        ss >> str2;
        str += str2;

        std::ofstream output;
        output.open(str.c_str(), std::ios::out);

        Output_structure structure_out;
        structure_out.output_from_cell_to_file(output, axis, type, 
            order, position);
        output << " axis_change = " << axis_change << std::endl;
        output << " n_type = ";
        for (int j = 0; j < n_type.size(); ++j){
            output << n_type[j] << " ";
        }
        output << std::endl;

        output.close();
    }

}

void Output::output_cluster_out(
    const vector2_pair& all_unique_cluster_output, 
    const vector2_i& all_unique_cluster, 
    const dmatrix& axis_primitive, const std::vector<int>& n_atoms, 
    const std::vector<dvector>& position_primitive, 
    const int& n_sublattice, const dmatrix& distance_array){

    std::ofstream output;
    output.open("cluster.out", std::ios::out);

    output.precision(15);
    output.width(17);
    output << " n_sublattice = " << n_sublattice << std::endl;

    Output_structure os;
    os.output_from_structure_to_file(output, axis_primitive,
            n_atoms, position_primitive);

    int n_count(1);
    output << " # index  nbody  positions " << std::endl;
    for (int i = 0; i < all_unique_cluster_output.size(); ++i){
        output << n_count << " " << all_unique_cluster_output[i].size() << " ";
        for (int j = 0; j < all_unique_cluster_output[i].size(); ++j){
            output << all_unique_cluster_output[i][j].first << " ";
            output << all_unique_cluster_output[i][j].second + 1 << " ";
        }
        if (all_unique_cluster_output[i].size() == 2){
            int atom1 = all_unique_cluster[i][0];
            int atom2 = all_unique_cluster[i][1];
            output << " # distance = " 
                << distance_array(atom1,atom2) << std::endl;
        }
        else{
            output << std::endl;
        }
        ++n_count;
    }

    output.close();

}

void Output::output_correlation(const dvector& correlation_array, 
    const std::vector<int>& n_clusters,
    const std::vector<dmatrix>& cluster_function_array, 
    const std::vector<int>& cluster_number_array, 
    const char* structure_name,
    const std::vector<int>& n_type){

    std::ofstream output;
    output.open("correlation.out", std::ios::out);

    output << int(correlation_array.size()) + 1
        << " # " << structure_name << std::endl;
    output << "0 0 1.000000000"<< std::endl;
    for (int i = 0; i < correlation_array.size(); ++i){
        output << i+1 << " " << cluster_number_array[i] + 1 << " ";
        double correlation = correlation_array(i);
        if (fabs(correlation) < 1e-12){
            correlation = 0;
        }
        output.precision(15);
        output.setf(std::ios::showpoint);
        output << correlation << " " << n_clusters[i] << std::endl;
    }
    output.close();

    output.open("cluster_function.out", std::ios::out);

    output << " n_type = ";
    for (int i = 0; i < n_type.size(); ++i){
        output << n_type[i] << " ";
    }
    output << std::endl;

    for (int i = 0; i < correlation_array.size(); ++i){
        output << i+1 << " " << cluster_number_array[i] + 1 << " ";
        output.setf(std::ios::showpoint);
        output.precision(15);
        output << cluster_function_array[i] << std::endl;
    }
    output.close();
}

void Output::output_lsf(const std::vector<int>& cluster_index,
    const dvector& eci, const double& square_error, 
    const double& cv_score, const dmatrix& inverse, 
    const int& inverse_check){

    std::ofstream output;
    output.open("eci", std::ios::out);

    output << " cluster_index = ";
    for (int i = 0; i < cluster_index.size(); ++i){
        output << cluster_index[i] << " ";
    }
    output << std::endl;
    output.precision(15);
    output.setf(std::ios::showpoint);
    output << " eci = " << eci << std::endl;
    output << " squared error = " <<  square_error << std::endl;
    output << " cv score = " << cv_score << std::endl;

    std::cout.precision(15);
    std::cout.setf(std::ios::showpoint);
    for (int i = 0; i < cluster_index.size(); ++i){
        std::cout << " eci " << cluster_index[i] 
            << "  " << eci(i) << std::endl;
    }
    std::cout << " squared error = " <<  square_error << std::endl;
    std::cout << " cv score = " << cv_score << std::endl;
    std::cout << " ### diagonal terms of precision matrix ###" << std::endl;
    if (inverse_check == 1){
        for (int i = 0; i < inverse.size1(); ++i){
            if (cluster_index[i] != 0){
                std::cout << " " << cluster_index[i] 
                    << "  " << inverse(i,i) << std::endl;
            }
        }
    }

    output.close();
}

void Output::output_cv_casp(const std::vector<double>& cv_casp, 
        const std::vector<int> n_structure_in_group, const double& cv_score){

    for (int i = 0; i < cv_casp.size(); ++i){
        std::cout << "  " << cv_casp[i] << " # group " << i+1 
            << " (n_structure = " << n_structure_in_group[i] << ")" 
            << std::endl;
    }

    std::cout << " total cv score = " << cv_score << std::endl;

}

void Output::output_ga(const vector2_i& pop_array, 
    const std::vector<double>& score_array, const int& iter){

    std::ofstream output;
    if (iter == -1){
        output.open("ga.out", std::ios::out);
    }
    else {
        output.open("ga.out", std::ios::app);
    }

    output << "  generation " << iter + 1 << std::endl;
    for (int i = 0; i < pop_array.size(); ++i){
        output << " cv_score = " << score_array[i] << "; cluster_index = ";
        for (int j = 0; j < pop_array[i].size(); ++j){
            output << pop_array[i][j] << " ";
        }
        output << std::endl;
    }

    output.close();

}

void Output::output_gss(const std::vector<double>& composition_array, 
    const std::vector<double>& energy_array,
    const std::vector<double>& error_array,
    const double& rms_error,
    const std::vector<double>& gs_comp,
    const std::vector<double>& gs_energy,
    const std::vector<int>& gs_line){

    std::ofstream output_gs;
    output_gs.precision(15);
    output_gs.setf(std::ios::showpoint);
    output_gs.open("gs.out", std::ios::out);
    output_gs << " # composition  energy  line(correlation) " << std::endl;
    for (int i = 0; i < gs_comp.size(); ++i){
        output_gs << gs_comp[i] << " " << gs_energy[i]
            << " " << gs_line[i] << std::endl;
    }

    output_gs.close();

    std::ofstream output;
    output.open("all_energy.out", std::ios::out);

    output.precision(7);
    output.setf(std::ios::showpoint);
    output << " # energy   composition " << std::endl;
    for (int i = 0; i < composition_array.size(); ++i){
        output << " " << energy_array[i] 
            << "  " << composition_array[i] << std::endl;
    }
    output.close();

    if (error_array.size() > 0){
        std::ofstream output2;
        output2.open("error.out", std::ios::out);

        output2.precision(7);
        output2.setf(std::ios::showpoint);
        output2 << " # rms of CE error = " << rms_error << std::endl;
        output2 << " # ce_energy   diff " << std::endl;
        for (int i = 0; i < energy_array.size(); ++i){
            output2 << " " << energy_array[i] 
                << "  " << error_array[i] << std::endl;
        }
        output2.close();
    }
    
}

void Output::output_mc
(std::vector<int>& spin, const double& energy_fin, const double& energy_ave,
 const double& specific_heat, const dvector& correlation_ave, 
 const std::vector<double>& temp_array, const dmatrix& axis, 
 const std::vector<int>& n_atoms, const std::vector<dvector>& position, 
 const std::vector<int>& cluster_index, const std::vector<int>& n_type, 
 std::vector<int>& spin_array){

    std::cout << "  energy of final structure = " << energy_fin << std::endl;
    std::cout << "  average energy = " << energy_ave << std::endl;
    std::cout << "  specific heat = " << specific_heat << std::endl;

    std::ofstream output;
    output.open("mc.out", std::ios::out);

    output << "  temperature = " << *(temp_array.end() - 1) << std::endl;
    output.precision(10);
    output << "  average energy = " << energy_ave << std::endl;
    output << "  specific heat = " << specific_heat << std::endl;
    output << "  average structure " << std::endl;
    output << "   c.f. of cluster " << cluster_index[0] << " = 1" << std::endl;
    for (int i = 0; i < correlation_ave.size(); ++i){
        output << "   c.f. of cluster " << cluster_index[i+1] << " = " 
            << correlation_ave(i) << std::endl;
    }
    output << "  energy of final structure = " << energy_fin << std::endl;
    output.close();

    std::ofstream output_structure;
    output_structure.open("structure_final", std::ios::out);

    mc_spin_to_structure(n_type, n_atoms, spin, spin_array);

    Output_structure structure_out;
    structure_out.output_from_cell_to_file
        (output_structure, axis, spin, spin_array, position);
    output_structure.close();

}

void Output::output_sqs
(const dvector& correlation_array, 
 const std::vector<int>& cluster_index, 
 const dmatrix& axis, const std::vector<int>& n_atoms, 
 const std::vector<dvector>& position, 
 const std::vector<int>& n_type, 
 std::vector<int>& spin, std::vector<int>& spin_array){

    std::ofstream output;
    output.open("sqs.out", std::ios::out);

    output.precision(10);
    output << "  Special quasirandom structure " << std::endl;
    for (int i = 0; i < correlation_array.size(); ++i){
        output << "   c.f. of cluster " << cluster_index[i] << " = " 
            << correlation_array(i) << std::endl;
    }
    output.close();

    std::ofstream output_structure;
    output_structure.open("structure_sqs", std::ios::out);

    mc_spin_to_structure(n_type, n_atoms, spin, spin_array);

    Output_structure structure_out;
    structure_out.output_from_cell_to_file
        (output_structure, axis, spin, spin_array, position);

    output_structure.close();

}

void Output::mc_spin_to_structure
(const std::vector<int>& n_type, const std::vector<int>& n_atoms,
 std::vector<int>& spin, std::vector<int>& spin_array){

    int n(0), n2(0);
    for (int i = 0; i < n_type.size(); ++i){
        int add = (i+1) * 1e5;
        for (int j = 0; j < n_type[i]; ++j){
            spin_array[n] += add;
            ++n;
        }
        for (int j = 0; j < n_atoms[i]; ++j){
            spin[n2] += add;
            ++n2;
        }
    }
    for (int i = n_type.size(); i < n_atoms.size(); ++i){
        spin_array.push_back(i+1000000);
        for (int j = 0; j < n_atoms[i]; ++j){
            spin.push_back(i+1000000);
        }
    }
}
