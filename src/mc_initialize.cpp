/******************************************************************************

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

	    Class for performing initialization for MC simulations
	
*****************************************************************************/

#include "mc_initialize.h"

MC_initialize::MC_initialize(const Input& ip){

    pi = 3.14159265358979323846;

    // construction of simulation supercell
    const dmatrix axis_primitive = ip.get_axis_primitive();
    const std::vector<dvector> position_primitive = ip.get_position_primitive();
    const std::vector<int> n_atoms_primitive = ip.get_n_atoms_primitive();
    const std::vector<int> type_primitive = ip.get_type_primitive();
    const imatrix supercell_matrix = ip.get_axis_change();

    Structure structure(axis_primitive, position_primitive, n_atoms_primitive);
    structure.supercell(supercell_matrix);

    axis = structure.get_axis();
    position = structure.get_position();
    n_atoms = structure.get_n_atoms();

    dmatrix inverse_supercell_matrix;
    math::invert(supercell_matrix, inverse_supercell_matrix);

    // number of steps
    const std::vector<int> n_type = ip.get_n_type();
    int n_atom_simulation 
        = accumulate(n_atoms.begin(), n_atoms.begin() + n_type.size(), 0);
    n_step_ann = ip.get_n_step_ann() * n_atom_simulation;
    n_step_eqv = ip.get_n_step_eqv() * n_atom_simulation;

    std::cout << "  n-atoms (simulation) = " << n_atom_simulation << std::endl;
    if (ip.get_simulation_type() == "sqs"){
        std::cout << "  n-steps = " << n_step_ann << std::endl;
    }
    else {
        std::cout << "  n-steps for annealing = " << n_step_ann << std::endl;
        std::cout << "  n-steps for averaging = " << n_step_eqv << std::endl;
    }

    // setting sublattice
    for (int i = 0; i < n_type.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            sublattice.push_back(i);
        }
    }
 
    // set eci, cluster_index, unique_cluster_array, cluster_function_array
    vector2_pair unique_cluster_array;
    sort_eci_cluster(ip, unique_cluster_array);
    eci *= supercell_matrix(0,0) * supercell_matrix(1,1) 
        * supercell_matrix(2,2);
    eci_without_empty *= supercell_matrix(0,0) 
        * supercell_matrix(1,1) * supercell_matrix(2,2);

    // preparing initial structure
    const std::vector<int> spin_array = ip.get_spin_array();
    const double symprec = ip.get_symprec();
    std::vector<double> complete_comp;
    if (ip.get_init_structure() == "none" or ip.get_simulation_type() == "sqs"){
        complete_comp = get_complete_comp_array(ip.get_comp_array(), n_type);
        spin = init_structure_rand
            (n_atoms, n_type, complete_comp, spin_array, symprec);
    }
    else {
        spin = init_structure_file
            (ip.get_init_structure(), position, supercell_matrix, 
             n_atom_simulation, spin_array, symprec);
    }

    // for calculation of electrostatic energy
    energy = 0.0;
    if (ip.get_i_ewald() == 1){
        // charge of initial structure
        charge = structure_to_charge
            (spin, spin_array, n_type, n_atoms, ip.get_charge_array());
        double total_charge = accumulate(charge.begin(), charge.end(), 0.0);
        if (fabs(total_charge) > 1e-5){
            std::cerr << " Total charge is not zero." << std::endl;
            exit(8);
        }
        ewald_obj = Ewald 
            (axis, position, structure.get_volume(), 
             structure.get_reciprocal_axis(), 
             charge, ip.get_eta(), ip.get_rmax(), ip.get_gmax());
        ewald_obj.calc_energy();
        energy = ewald_obj.get_energy() / ip.get_epsilon();
    }

    // calculation of permutation for calculating correlation functions, 
    // number of clusters and coordination around atoms
    Sym sym(axis_primitive, position_primitive, type_primitive, symprec);
    sym.supercell_symmetry(supercell_matrix);

    Correlation correlation;
    std::vector<imatrix> permutation_array = correlation.get_permutation_array
        (unique_cluster_array, position_primitive, position, 
         inverse_supercell_matrix, sym.get_double_rotate_matrix(),
         sym.get_trans_vector(), symprec);

    n_cluster_array = dvector(permutation_array.size());
    for (int i = 0; i < permutation_array.size(); ++i){
        n_cluster_array(i) = double(permutation_array[i].size1());
    }

    // setting coordination around atoms
    int binary(0);
    if (n_type.size() == count(n_type.begin(), n_type.end(), 2)) binary = 1;
    if (binary == 1) set_coordination_binary(permutation_array);
    else set_coordination(permutation_array);

    // correlation functions and energy for initial structure
    correlation_array = correlation.get_correlation_function
        (permutation_array, cluster_function_array, spin);
    permutation_array.clear();

    if (ip.get_simulation_type() == "cmc" 
            or ip.get_simulation_type() == "gcmc"){
        energy += eci(0) 
            + ublas::inner_prod(eci_without_empty, correlation_array);
    }
    else if (ip.get_simulation_type() == "sqs"){
        disorder_cf = calc_disorder_cf
            (complete_comp, unique_cluster_array, 
             n_atoms_primitive, n_type, spin_array);
    }
}

MC_initialize::~MC_initialize(){}

// set eci, cluster_index, unique_cluster_array, cluster_function_array
void MC_initialize::sort_eci_cluster
(const Input& ip, vector2_pair& unique_cluster_array){

    const std::string simulation_type = ip.get_simulation_type();
    std::vector<int> cluster_number_array = ip.get_cluster_number_array();
    cluster_index = ip.get_cluster_index();
    unique_cluster_array = ip.get_unique_cluster_array();
    cluster_function_array = ip.get_cluster_function_array();

    if (simulation_type == "cmc" or simulation_type == "gcmc"){
        eci = ip.get_eci();
        eci_without_empty = dvector(eci.size()-1);
        dvector eci_tmp(eci.size());
        std::vector<int> cluster_index_tmp(cluster_index.size());
        int n = 1;
        for (int i = 0; i < cluster_index.size(); ++i){
            if (cluster_index[i] > 0){
                eci_tmp(n) = eci(i);
                eci_without_empty(n-1) = eci(i);
                cluster_index_tmp[n] = cluster_index[i];
                ++n;
            }
            else {
                eci_tmp(0) = eci(i);
                cluster_index_tmp[0] = cluster_index[i];
            }
        }
        eci = eci_tmp;
        cluster_index = cluster_index_tmp;
    }

    vector2_pair unique_cluster_array_tmp;
    std::vector<dmatrix> cluster_function_array_tmp;
    for (int i = 0; i < cluster_index.size(); ++i){
        int serial_number = cluster_index[i];
        if (serial_number > 0){
            int cluster_number = cluster_number_array[serial_number-1];
            unique_cluster_array_tmp.push_back
                (unique_cluster_array[cluster_number-1]);
            cluster_function_array_tmp.push_back
                (cluster_function_array[serial_number-1]);
        }
    }
    unique_cluster_array = unique_cluster_array_tmp;
    cluster_function_array = cluster_function_array_tmp;
}

std::vector<int> MC_initialize::init_structure_rand
(const std::vector<int>& n_atoms, const std::vector<int>& n_type, 
 const std::vector<double>& complete_comp_array, 
 const std::vector<int>& spin_array, const double& symprec){

    std::vector<int> spin_structure;
    srand(clock()+getpid());

    int count(0); double n_atom;
    for (int i = 0; i < n_type.size(); ++i){
        std::vector<int> spin_tmp;
        for (int j = 0; j < n_type[i]; ++j){
            n_atom = n_atoms[i] * complete_comp_array[count];
            if (fabs(n_atom - round(n_atom)) > symprec){
                std::cerr << " error : n-atoms * composition is "
                    << "not an integer. " << std::endl;
                exit(8);
            }
            for (int k = 0; k < round(n_atom); ++k){
                spin_tmp.push_back(spin_array[count]);
            }
            ++count;
        }
        while (spin_tmp.size() != 0){
            int n_rand = rand() % spin_tmp.size();
            spin_structure.push_back(spin_tmp[n_rand]);
            spin_tmp.erase(spin_tmp.begin()+n_rand);
        }
    }
    return spin_structure;
}

std::vector<int> MC_initialize::init_structure_file
(const std::string& init_structure, const std::vector<dvector>& position, 
 const imatrix& supercell_matrix, const int& n_atom_simulation, 
 const std::vector<int>& spin_array, const double& symprec){

    std::vector<int> spin_structure;

    Parse_structure structure(init_structure.c_str());
    dmatrix axis_structure = structure.get_axis();
    std::vector<int> n_atoms_structure = structure.get_num_atoms();
    std::vector<dvector> position_structure = structure.get_position();

    dmatrix axis_change_structure;
    std::vector<int> spin, n_type_structure;

    Parse_input parse_structure(init_structure.c_str());
    parse_structure.assignNeed("axis_change", axis_change_structure);
    parse_structure.assignNeed("n_type", n_type_structure);

    int n_type_sum 
        = accumulate(n_type_structure.begin(), n_type_structure.end(), 0);
    for (int i = 0; i < n_type_sum; ++i){
        for (int j = 0; j < n_atoms_structure[i]; ++j){
            spin.push_back(spin_array[i]);
        }
    }

    dmatrix inverse_axis_change_structure;
    math::invert(axis_change_structure, inverse_axis_change_structure);
    dmatrix change = ublas::prod
        (inverse_axis_change_structure, supercell_matrix);

    for (int i = 0; i < n_atom_simulation; ++i){
        dvector pos = ublas::prod(change, position[i]);
        for (int j = 0; j < 3; ++j){
            while (pos(j) > 1 - symprec){
                pos(j) -= 1.0;
            }
            while (pos(j) < -symprec){
                pos(j) += 1.0;
            }
        }
        for (int j = 0; j < position_structure.size(); ++j){
            if (ublas::norm_inf(pos - position_structure[j]) < symprec){
                spin_structure.push_back(spin[j]);
                break;
            }
        }
    }
    return spin_structure;
}

void MC_initialize::set_coordination_binary
(std::vector<imatrix>& permutation_array){

    // coordination
    int n_atom_simulation = spin.size();
    coordination.resize(n_atom_simulation);
    n_cluster_coordination.resize(n_atom_simulation);
    for (int i = 0; i < n_atom_simulation; ++i){
        coordination[i].resize(permutation_array.size());
        n_cluster_coordination[i].resize(permutation_array.size());
    }

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int i = 0; i < permutation_array.size(); ++i){
        std::multiset<std::vector<int> > cluster_array;
        for (int j = 0; j < permutation_array[i].size1(); ++j){
            std::multiset<int> cluster_set;
            for (int k = 0; k < permutation_array[i].size2(); ++k){
                cluster_set.insert(permutation_array[i](j,k));
            }
            std::vector<int> cluster(cluster_set.size());
            std::copy(cluster_set.begin(), cluster_set.end(), cluster.begin());
            cluster_array.insert(cluster);
        }
        std::multiset<std::vector<int> >::iterator it = cluster_array.begin();
        while (it != cluster_array.end()){
            std::vector<int> cluster = *it;
            int n_same = cluster_array.count(*it);
            std::set<int> cluster_duplication_check;
            for (int j = 0; j < cluster.size(); ++j){
                cluster_duplication_check.insert(cluster[j]);
            }
            std::set<int>::iterator it_set = cluster_duplication_check.begin();
            while (it_set != cluster_duplication_check.end()){
                coordination[*it_set][i].push_back(cluster);
                n_cluster_coordination[*it_set][i].push_back(n_same);
                ++it_set;
            }
            for (int j = 0; j < n_same; ++j){
                ++it;
            }
        }
    }
}

void MC_initialize::set_coordination
(std::vector<imatrix>& permutation_array){

    // coordination
    int n_atom_simulation = spin.size();
    coordination.resize(n_atom_simulation);
    for (int i = 0; i < n_atom_simulation; ++i){
        coordination[i].resize(permutation_array.size());
    }

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int i = 0; i < permutation_array.size(); ++i){
        for (int j = 0; j < permutation_array[i].size1(); ++j){
            std::vector<int> cluster;
            for (int k = 0; k < permutation_array[i].size2(); ++k){
                cluster.push_back(permutation_array[i](j,k));
            }
            for (int k = 0; k < cluster.size(); ++k){
                coordination[cluster[k]][i].push_back(cluster);
            }
        }
    }
}

dvector MC_initialize::calc_disorder_cf
(const std::vector<double>& complete_comp_array, 
 const vector2_pair& unique_cluster_array,
 const std::vector<int>& n_atoms_primitive,
 const std::vector<int>& n_type,
 const std::vector<int>& spin_array){

    dvector output(correlation_array.size());
    vector2_d average_spin;
    int n(0);
    for (int i = 0; i < n_type.size(); ++i){
        std::vector<double> average_array;
        int start = accumulate(n_type.begin(), n_type.begin() + i, 0);
        for (int j = 0; j < n_type[i]; ++j){
            double average(0.0);
            for (int k = 0; k < n_type[i]; ++k){
                average += pow(spin_array[start+k], j) 
                    * complete_comp_array[start+k];
            }
            average_array.push_back(average);
        }
        average_spin.push_back(average_array);
    }

    vector2_i cluster_sublattice;
    for (int i = 0; i < unique_cluster_array.size(); ++i){
        std::vector<int> sublattice;
        for (int j = 0; j < unique_cluster_array[i].size(); ++j){
            int atom = unique_cluster_array[i][j].second;
            for (int k = 0; k < n_atoms_primitive.size(); ++k){
                atom -= n_atoms_primitive[k];
                if (atom < 0){
                    sublattice.push_back(k);
                    break;
                }
            }
        }
        cluster_sublattice.push_back(sublattice);
    }

    for (int i = 0; i < cluster_function_array.size(); ++i){
        dmatrix function = cluster_function_array[i];
        if (function.size1() < 2){
            std::cout << " Error: Empty or point clusters are "
                << "included in cluster_index." << std::endl;
            exit(8);
        }
        double prod(1.0);
        for (int j = 0; j < function.size1(); ++j){
            double sum(0.0);
            int sub = cluster_sublattice[i][j];
            for (int k = 0; k < function.size2(); ++k){
                sum += function(j,k) * average_spin[sub][k];
            }
            prod *= sum;
        }
        output(i) = prod;
    }
    return output; 
}

std::vector<double> MC_initialize::get_complete_comp_array
(const std::vector<double>& comp_array, const std::vector<int>& n_type){

    std::vector<double> complete_comp_array;

    int n(0), n_type_sum(0);         
    for (int i = 0; i < n_type.size(); ++i){
        double comp_first(1.0);
        for (int j = 0; j < n_type[i] - 1; ++j){
            complete_comp_array.push_back(comp_array[n]);
            comp_first -= comp_array[n];
            ++n;
        }         
        complete_comp_array.insert
            (complete_comp_array.begin() + n_type_sum, comp_first);
        n_type_sum += n_type[i];
    }

    return complete_comp_array;
}

std::vector<double> MC_initialize::structure_to_charge
(const std::vector<int>& spin, const std::vector<int>& spin_array, 
 const std::vector<int>& n_type, const std::vector<int>& n_atoms,
 const std::vector<double>& charge_array){

    std::vector<double> output;

    int n(0), start;
    for (int i = 0; i < n_type.size(); ++i){
        start = accumulate(n_type.begin(), n_type.begin() + i, 0);
        for (int j = 0; j < n_atoms[i]; ++j){
            for (int k = start; k < start + n_type[i]; ++k){
                if (spin[n] == spin_array[k]){
                    output.push_back(charge_array[k]);
                }
            }
            ++n;
        }
    }
    start = accumulate(n_type.begin(), n_type.begin() + n_type.size(), 0);
    for (int i = n_type.size(); i < n_atoms.size(); ++i){
        int charge_index = start + i - n_type.size();
        for (int j = 0; j < n_atoms[i]; ++j){
            output.push_back(charge_array[charge_index]);
        }
    }

    return output;
}

// structure
const dmatrix& MC_initialize::get_axis() const{
    return axis;
}
const std::vector<int>& MC_initialize::get_n_atoms() const{
    return n_atoms;
}
const std::vector<dvector>& MC_initialize::get_position() const{
    return position;
}

// size = number of atoms
const std::vector<int>& MC_initialize::get_spin() const{
    return spin;
}
const std::vector<int>& MC_initialize::get_sublattice() const{
    return sublattice;
}

// input information
const int& MC_initialize::get_n_step_ann() const{
    return n_step_ann;
}
const int& MC_initialize::get_n_step_eqv() const{
    return n_step_eqv;
}

// for sqs exploration
const dvector& MC_initialize::get_disorder_cf() const{
    return disorder_cf;
}

// for c.f. calculation
const dvector& MC_initialize::get_eci_without_empty() const{
    return eci_without_empty;
}
const std::vector<int>& MC_initialize::get_cluster_index() const{
    return cluster_index; 
}
const std::vector<dmatrix>& MC_initialize::get_cluster_function_array() const{
    return cluster_function_array;
}
const dvector& MC_initialize::get_n_cluster_array() const{
    return n_cluster_array;
}
const vector4_i& MC_initialize::get_coordination() const{
    return coordination;
}
const vector3_i& MC_initialize::get_n_cluster_coordination() const{
    return n_cluster_coordination;
}

// energy and c.f. of initial structure
const double& MC_initialize::get_energy() const{
    return energy;
}
const dvector& MC_initialize::get_correlation_array() const{
    return correlation_array;
}

const std::vector<double>& MC_initialize::get_charge() const{
    return charge;
}
const Ewald& MC_initialize::get_ewald_obj() const{
    return ewald_obj;
}

