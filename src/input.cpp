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

	    Class for reading input files
	
*****************************************************************************/

#include "input.h"

Input::Input(){

    pi = 3.14159265358979323846;
}

Input::~Input(){}

void Input::set_position_simulation
(const int& n_sublattice,
 const std::vector<dvector>& position_all, 
 const std::vector<int>& n_atoms){

    int n(0);
    for (int i = 0; i < n_sublattice; ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            position_simulation.push_back(position_all[n]);
            ++n;
        }
    }
}

void Input::input_derivative(){

    Parse_input parse("derivative.in");

    parse.assignNeed("n_lattice", n_lattice);
    parse.assignNeed("n_type", n_type);
    parse.assign("symprec", symprec, 1e-5);
    n_sublattice = n_type.size();

    Parse_structure punit("primitive_lattice");
    axis_primitive = punit.get_axis();
    n_atoms_primitive = punit.get_num_atoms();
    position_primitive = punit.get_position();
    type_primitive = punit.get_type();

    for (int i = 0; i < n_atoms_primitive.size(); ++i){
        n_atoms.push_back(n_atoms_primitive[i] * n_lattice);
        if (n_atoms[i] > 32 and i < n_type.size()){
            std::cout << " Warning: n_atoms[" << i << "],"
                << " derivative may not work." << std::endl;
        }
    }

    // for fixed composition
    std::vector<double> comp_array, comp_array_default;
    parse.assign("comp", comp_array, comp_array_default);
    if (comp_array.size() > 0){
        int comp_size(0);
        for (int i = 0; i < n_type.size(); ++i){
            comp_size += n_type[i] - 1;
        }
        if (comp_size == comp_array.size()){
            int comp_index(0);
            for (int i = 0; i < n_type.size(); ++i){
                double comp_first(1.0);
                std::vector<double> tmp;
                for (int j = 0; j < n_type[i] - 1; ++j){
                    comp_first -= comp_array[comp_index];
                    tmp.push_back(comp_array[comp_index]);
                    ++comp_index;
                }
                std::vector<double> comp_all;
                comp_all.push_back(comp_first);
                comp_all.insert(comp_all.end(), tmp.begin(), tmp.end());
                for (int j = 0; j < comp_all.size(); ++j){
                    double n = comp_all[j] * n_atoms[i];
                    if (fabs(n - round(n)) < symprec){
                        n_atoms_derivative.push_back(round(n));
                    }
                    else {
                        std::cerr << " n_atoms (n_atoms of prim. * n_lattice) "
                        "* composition != integer" << std::endl;
                        exit(8);
                    }
                }
            }
        }
        else {
            std::cerr << " number of elements in comp should be " 
                << comp_size << " or 0 for searching derivative "
                << "structure in a whole range of compositions. "<< std::endl;
            exit(8);
        }
    }

    set_position_simulation(n_sublattice, position_primitive, 
        n_atoms_primitive);

    Error_check error;
    error.check_derivative(n_atoms_primitive, n_type);

    // standard output
    std::cout << "  derivative structures with " << n_lattice 
        << " primitive lattices " << std::endl;
    int n(0), index(0);
    for (int i = 0; i < n_type.size(); ++i){
        std::cout << "  sublattice " << i+1 << " : n_type = " 
            << n_type[i] << std::endl;
        if (n_atoms_derivative.size() > 0){
            std::cout << "  ";
            for (int j = 0; j < n_type[i]; ++j){
                std::cout << " type " << index << " = "
                << n_atoms_derivative[index] << " atoms; ";
                ++index;
            }
            std::cout << std::endl;
        }

        for (int j = 0; j < n_atoms_primitive[i]; ++j){
            std::cout << "   position" << n+1 << " = " 
                << position_primitive[n] << std::endl;
            ++n;
        }

    }
}

void Input::input_derivative_to_structure(){

    Parse_input parse("derivative.out");
    parse.assignNeed("n_type", n_type);

    Parse_structure punit("derivative.out");
    axis_primitive = punit.get_axis();
    n_atoms_primitive = punit.get_num_atoms();
    position_primitive = punit.get_position();
    type_primitive = punit.get_type();

    int n_all_atoms(0);
    for (int i = 0; i < n_atoms_primitive.size(); ++i){
        n_all_atoms += n_atoms_primitive[i];
    }

    std::ifstream input("derivative.out");
    std::string line;
    for (int i = 0; i < n_all_atoms + 6; ++i){
        std::getline(input, line);
    }

    while (std::getline(input, line) && !input.eof()){
        std::stringstream ss;
        int index;
        std::string label;
        imatrix hnf, snf, left;
        ss << line;
        ss >> index >> label >> hnf >> snf >> left;
        hnf_array.push_back(hnf);
        snf_array.push_back(snf);
        left_array.push_back(left);

        std::vector<int> tmp;
        for (int i = 0; i < label.size(); ++i){
            int atom_type;
            std::stringstream ss2;
            ss2 << label[i];
            ss2 >> atom_type;
            tmp.push_back(atom_type);
        }
        label_array.push_back(tmp);
    }
}

void Input::input_derivative_to_structure_fix_comp(){

    Parse_input parse("derivative.out");
    parse.assignNeed("n_type", n_type);

    Parse_structure punit("derivative.out");
    axis_primitive = punit.get_axis();
    n_atoms_primitive = punit.get_num_atoms();
    position_primitive = punit.get_position();
    type_primitive = punit.get_type();

    int n_all_atoms(0);
    for (int i = 0; i < n_atoms_primitive.size(); ++i){
        n_all_atoms += n_atoms_primitive[i];
    }

    std::ifstream input("derivative.out");
    std::string line;
    for (int i = 0; i < n_all_atoms + 6; ++i){
        std::getline(input, line);
    }

    while (std::getline(input, line) && !input.eof()){
        std::stringstream ss;
        int index;
        std::string label;
        imatrix hnf, snf, left;
        ss << line;
        ss >> index >> label >> hnf >> snf >> left;

        std::vector<int> tmp;
        for (int i = 0; i < label.size(); ++i){
            int atom_type;
            std::stringstream ss2;
            ss2 << label[i];
            ss2 >> atom_type;
            tmp.push_back(atom_type);
        }

        int n_count = static_cast<int>(count(tmp.begin(), tmp.end(), 1));
        if (n_count == 2){
            hnf_array.push_back(hnf);
            snf_array.push_back(snf);
            left_array.push_back(left);
            label_array.push_back(tmp);
        }
    }
}

void Input::input_cluster(){

    Parse_input parse("cluster.in");
    parse.assign("max_cluster", max_cluster, 4);
    parse.assignNeed("trunc_distance", trunc_distance);
    parse.assign("symprec", symprec, 1e-5);
    parse.assign("n_sublattice", n_sublattice, 1);

    Parse_structure punit("primitive_lattice");
    axis_primitive = punit.get_axis();
    n_atoms_primitive = punit.get_num_atoms();
    position_primitive = punit.get_position();
    type_primitive = punit.get_type();

    set_position_simulation(n_sublattice, position_primitive, 
        n_atoms_primitive);
}

void Input::input_correlation(int argc, char *argv[]){

    read_cluster_out();

    std::string file = read_file_foption(argc, argv);

    Parse_structure structure(file.c_str());
    axis = structure.get_axis();
    n_atoms = structure.get_num_atoms();
    position = structure.get_position();
    type = structure.get_type();

    Parse_input parse(file_name);
    parse.assignNeed("axis_change", axis_change);
    parse.assignNeed("n_type", n_type);

    Parse_input parse_correlation("correlation.in");
    parse_correlation.assignNeed("spin", spin_array);
    parse_correlation.assign("symprec", symprec, 1e-5);

    Error_check error;
    error.check_correlation(axis_primitive, axis_change, axis, 
        n_atoms_primitive, n_atoms, n_type, spin_array, symprec);

    std::vector<int> type_tmp;
    int n(0);
    for (int i = 0; i < n_type.size(); ++i){
        for (int j = 0; j < n_type[i]; ++j){
            for (int k = 0; k < n_atoms[n]; ++k){
                spin_structure.push_back(spin_array[n]);
                type_tmp.push_back(i);
            }
            ++n;
        }
    }

    int n_type_sum = accumulate(n_type.begin(), n_type.end(), 0);
    for (int i = n_type_sum; i < n_atoms.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            spin_structure.push_back(1e4);
            type_tmp.push_back(i);
        }
    }
    type = type_tmp;

    // standard output
    for (int i = 0; i < n_type_sum; ++i){
        std::cout << "  spin value of " << i + 1 << "th atom type is " 
            << spin_array[i] << std::endl;
    }
}

void Input::input_lsf(){

    read_energy("dft_energy");
    read_correlation("correlation");

    Parse_input parse("lsf.in");
    set_weight(parse);
    parse.assignNeed("cluster_index", cluster_index);

    Error_check check;
    check.check_lsf_files(energy, correlation);

    dmatrix correlation_calc(energy.size(), cluster_index.size());
    for (int i = 0; i < energy.size(); ++i){
        for (int j = 0; j < cluster_index.size(); ++j){
            correlation_calc(i,j) = correlation(i,cluster_index[j]);
        }
    }
    correlation = correlation_calc;
}

void Input::input_cv_casp(int argc, char *argv[]){

    file_name = NULL;
    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"g:"))!=-1){
        switch( option ){
            case 'g' : file_name = optarg;
                       break;
        }
    }
    if (file_name == NULL){
        std::cerr << " File name of group is needed (option -g)." << std::endl;
        exit(8);
    }
    read_group(file_name);

    input_lsf();
}

void Input::input_ga(){

    read_energy("dft_energy");
    read_correlation("correlation");
    Error_check check;
    check.check_lsf_files(energy, correlation);

    Parse_input parse("ga.in");
    set_weight(parse);

    parse.assign("n_loop", nloop, 0);
    if (nloop > 0) get_temp_array(parse, 1);

    int min_index, max_index;
    parse.assign("min_index", min_index, 0);
    parse.assign("max_index", max_index, int(correlation.size2()-1));
    parse.assignNeed("n_cluster", n_cluster);

    std::vector<int> ex_array, empty_array;
    parse.assign("base_cluster", base_array, empty_array);
    parse.assign("ex_cluster", ex_array, empty_array);

    parse.assign("n_pop", n_pop, 30);
    parse.assign("max_iteration", max_iteration, 300);
    parse.assign("p_mating", p_mating, 0.9);
    parse.assign("p_mutation", p_mutation, 0.03);
    parse.assign("n_elite", n_elite, 1);

    for (int i = min_index; i < max_index + 1; ++i){
        if (find(ex_array.begin(), ex_array.end(),i) == ex_array.end()){
            cluster_candidate.push_back(i);
        }
    }
}

void Input::input_gss(){

    read_correlation("correlation");
    read_eci("eci");
    
    std::ifstream input("gss.in");
    if (input.fail()){
        composition_cluster = 1;
        composition_definition.push_back(0.5);
        composition_definition.push_back(-0.5);
    }
    else {
        Parse_input input("gss.in");
        input.assignNeed("comp_def", composition_definition);
        input.assignNeed("comp_cluster", composition_cluster);
        input.assign("error_check", error_check, 0);
    }

    if (error_check == 1){
        read_energy("dft_energy");
        Error_check check;
        check.check_lsf_files(energy, correlation);
    }
}

void Input::input_mc(){

    read_eci("eci"); // eci, cluster_index
    read_cluster_out(); // primitive lattice, n_sublattice, cluster positions
    read_cluster_function_out(); // n_type, cluster_function

    Parse_input input("mc.in");
    input.assign("symprec", symprec, 1e-5);

    // cell size, number of atoms and volume of simulation supercell 
    // (n_atom and volume is used for calculating electrostatic energy)
    Structure st(axis_primitive, position_primitive, n_atoms_primitive);
    volume = st.get_volume();
    int n_atom = position_primitive.size();

    std::vector<int> cell_expansion;
    input.assignNeed("cell_expansion", cell_expansion);
    axis_change = ublas::identity_matrix<double>(3);
    for (int i = 0; i < 3; ++i){
        axis_change(i,i) = cell_expansion[i];
        n_atom *= cell_expansion[i];
        volume *= cell_expansion[i];
    }

    // temperature and number of steps
    int i_temp;
    input.assign("i_temp", i_temp, 0);
    get_temp_array(input, i_temp);
    input.assign("n_step_ann", n_step_ann, 0);
    input.assignNeed("n_step_eqv", n_step_eqv);

    // initial structure
    input.assignNeed("spin", spin_array);
    error_spin_array_size(spin_array, n_type);

    const std::string init_structure_default("none"), 
        simulation_type_default("cmc");
    input.assign("init_structure", init_structure, init_structure_default);
    input.assign("simulation_type", simulation_type, simulation_type_default);

    if (init_structure == "none"){
        if (simulation_type == "cmc"){
            // comp(B,C,...); comp(A) = 1-comp(B)-comp(C)-...
            input.assignNeed("comp", comp_array); 
            error_comp_array_size(comp_array, spin_array, n_type);
        }
        else if (simulation_type == "gcmc"){
            for (int i = 0; i < n_type.size(); ++i){
                for (int j = 0; j < n_type[i] - 1; ++i){
                    comp_array.push_back(0.0);
                }
            }
        }
    }

    if (simulation_type == "gcmc"){
        // mu = mu(B) - mu(A); mu(C) - mu(A) ...
        // grand pot. = E - TS - mu * comp(B)
        input.assignNeed("mu", mu); 
        mu.insert(mu.begin(), 0.0);
        int pseudo;
        input.assign("pseudo", pseudo, 0); 
        if (pseudo > 1){
            for (int i = 0; i < pseudo; ++i){
                std::vector<int> array;
                std::string name("n_change"), index;
                std::stringstream ss;
                ss << i+1;
                ss >> index;
                input.assignNeed((name+index).c_str(), array); 
                n_change_array.push_back(array);
            }
        }
    }

    // for calculating electrostatic energy
    input.assign("i_ewald", i_ewald, 0);

    if (i_ewald == 1){
        read_ewald_input(input, n_atom, volume);
        input.assign("epsilon", epsilon, 1.0);
        input.assignNeed("charge", charge_array);
        if (charge_array.size() != 
                spin_array.size() + n_atoms_primitive.size() - n_type.size()){
            std::cerr << " Number of charge must be equal"
                << " to number of atom types." << std::endl;
            exit(8);
        }
    }
    else if (i_ewald != 0){
        std::cerr << " i_ewald must be 0 or 1. " << std::endl;
        exit(8);
    }
}

void Input::input_sqs(){

    simulation_type = "sqs";
    init_structure = "none";
    i_ewald = 0;

    read_cluster_out(); // primitive lattice, n_sublattice, cluster positions
    read_cluster_function_out(); // n_type, cluster_function

    Parse_input input("sqs.in");
    input.assign("symprec", symprec, 1e-5);
    input.assignNeed("cluster_index", cluster_index);

    // cell size of simulation supercell
    std::vector<int> cell_expansion;
    input.assignNeed("cell_expansion", cell_expansion);
    axis_change = ublas::identity_matrix<double>(3);
    for (int i = 0; i < 3; ++i){
        axis_change(i,i) = cell_expansion[i];
    }

    // temperature and number of steps
    get_temp_array(input, 1);
    input.assign("n_step", n_step_ann, 0);

    // initial structure
    input.assignNeed("spin", spin_array);
    error_spin_array_size(spin_array, n_type);

    // comp(B,C,...); comp(A) = 1-comp(B)-comp(C)-...
    input.assignNeed("comp", comp_array); 
    error_comp_array_size(comp_array, spin_array, n_type);

    std::string criterion_default("rms");
    input.assign("criterion", criterion, criterion_default);
    if (criterion != "rms" and criterion != "abs"){
        std::cerr << " criterion = rms or abs." << std::endl;
        exit(8);
    }
}

void Input::error_spin_array_size
(const std::vector<int>& spin_array, const std::vector<int>& n_type){

    if (spin_array.size() != accumulate(n_type.begin(), n_type.end(), 0)){
        std::cerr << " Sum of n_type must be equal to number of spins." 
            << std::endl;
        exit(8);
    }
}

void Input::error_comp_array_size
(const std::vector<double>& comp_array, 
 const std::vector<int>& spin_array, const std::vector<int>& n_type){

    if (spin_array.size() - n_type.size() != comp_array.size()){
        std::cerr << " Number of comp. must be " 
            << spin_array.size() - n_type.size() << "." << std::endl;
        exit(8);
    }
}

void Input::input_ewald(int argc, char *argv[]){

    std::string file = read_file_foption(argc, argv);

    Parse_structure structure(file.c_str());
    axis = structure.get_axis();
    n_atoms = structure.get_num_atoms();
    position = structure.get_position();
    volume = structure.get_volume();
    reciprocal_axis = structure.get_reciprocal_axis();

    Parse_input input("ewald.in");
    input.assignNeed("charge", charge_array);
    if (charge_array.size() != n_atoms.size()){
        std::cerr << " error: number of charge types "
            "is not equal to number of atom types. " << std::endl;
        exit(8);
    }
    for (int i = 0; i < n_atoms.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            charge.push_back(charge_array[i]);
        }
    }
    read_ewald_input(input, position.size(), volume);
}

void Input::input_ti(int argc, char *argv[]){

    file_name = NULL;
    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"f:"))!=-1){
        switch( option ){
            case 'f' : file_name = optarg;
                       break;
        }
    }

    if (file_name == NULL){
        std::cerr 
            << " Specify file name of structure (option -f)." << std::endl;
        exit(8);
    }

    std::cerr << "  caution: ti is available only for binary system. " 
        << std::endl;
    
    Parse_input parse("ti.in");
    parse.assign("initial_potential", init_potential, 0.0);

    std::ifstream input(file_name);
    if (input.fail()){
        std::cerr << "Error: Could not open "<< file_name 
            << " file." << std::endl;
        exit(8);
    }

    std::string line;
    while( std::getline(input, line) && !input.eof()){
        int n_data;
        std::stringstream str;
        str << line;
        str >> n_data;
        std::vector<double> array1, array2, array3, array4;
        for (int i = 0; i < n_data; ++i){
            std::stringstream str2;
            std::getline(input, line);
            str2 << line;
            double mu, temp, comp, energy;
            str2 >> mu >> temp >> comp >> energy;
            array1.push_back(mu);
            array2.push_back(temp);
            array3.push_back(comp);
            array4.push_back(energy);
        }
        ti_mu_array.push_back(array1);
        ti_temp_array.push_back(array2);
        ti_comp_array.push_back(array3);
        ti_energy_array.push_back(array4);
    }
}

void Input::read_cluster_out(){

    Parse_structure prim("cluster.out");
    axis_primitive = prim.get_axis();
    n_atoms_primitive = prim.get_num_atoms();
    position_primitive = prim.get_position();
    type_primitive = prim.get_type();

    int n_all_atoms_primitive(0);
    for (int i = 0; i < n_atoms_primitive.size(); ++i){
        n_all_atoms_primitive += n_atoms_primitive[i];
    }

    Parse_input parse("cluster.out");
    parse.assignNeed("n_sublattice", n_sublattice);

    std::ifstream input("cluster.out");
    std::string line;
    for (int i = 0; i < n_all_atoms_primitive + 6; ++i){
        std::getline(input, line);
    }

    while (std::getline(input, line) && !input.eof()){
        std::stringstream ss;
        ss << line;
        int cluster_number, n_cluster;
        ss >> cluster_number >> n_cluster;
        ivector lattice;
        int atom_index;
        std::vector<std::pair<ivector, int> > cluster_information;
        for (int i = 0; i < n_cluster; ++i){
            ss >> lattice >> atom_index;
            std::pair<ivector, int> tmp_pair;
            tmp_pair.first = lattice;
            tmp_pair.second = atom_index - 1;
            cluster_information.push_back(tmp_pair);
        }
        unique_cluster_array.push_back(cluster_information);
    }
}

void Input::read_group(const char* file_name){

    std::ifstream input(file_name);
    if (input.fail()){
        std::cerr << "Error: Could not open " << file_name << std::endl;
        exit(8);
    }

    std::string line;
    while (std::getline(input, line) && !input.eof()){
        int group_index;
        std::string structure_index;
        std::stringstream str;
        str << line;
        str >> structure_index >> group_index;
        structure_index_array.push_back(structure_index);
        group_index_array.push_back(group_index);
    }
}

void Input::read_energy(const char* file_name){

    std::ifstream input(file_name);
    if (input.fail()){
        std::cerr << "Error: Could not open " << file_name << std::endl;
        exit(8);
    }
    std::string line;
    std::vector<double> energy_array;
    while ( std::getline(input, line) && !input.eof()){
        double energy_tmp;
        std::string st_index;
        std::stringstream str;
        str << line;
        str >> energy_tmp >> st_index;
        line.erase();
        energy_array.push_back(energy_tmp);
        energy_structure_index.push_back(st_index);
    }

    energy.resize(energy_array.size());
    for (int i = 0; i < energy_array.size(); ++i){
        energy(i) = energy_array[i];
    }
}

void Input::read_correlation(const char* file_name){

    std::ifstream input(file_name);
    if (input.fail()){
        std::cerr << "Error: Could not open " << file_name << std::endl;
        exit(8);
    }

    vector2_d correlation_array;
    std::vector<double> correlation_vector;
    std::string line;
    while (std::getline(input,line) && !input.eof()){
        int n_all_cluster;
        std::stringstream str;
        str << line;
        str >> n_all_cluster;
        for (int i = 0; i < n_all_cluster; ++i){
            std::string line2;
            std::stringstream str2;
            double correlation;
            int cluster_serial_number, cluster_number;
            std::getline(input, line2);
            str2 << line2;
            str2 >> cluster_serial_number >> cluster_number >> correlation;
            correlation_vector.push_back(correlation);
        }
        correlation_array.push_back(correlation_vector);
        correlation_vector.clear();
    }

    correlation = dmatrix(correlation_array.size(), 
        correlation_array[0].size());
    for (int i = 0; i < correlation_array.size(); ++i){
        for (int j = 0; j < correlation_array[0].size(); ++j){
            correlation(i,j) = correlation_array[i][j];
        }
    }
}

void Input::read_weight(const char* file_name){

    std::ifstream input(file_name);
    if (input.fail()){
        std::cerr << "Error: Could not open " << file_name << std::endl;
        exit(8);
    }

    std::vector<double> weight_array;
    std::string line;
    while ( std::getline(input, line) && !input.eof()){
        double weight_tmp;
        std::stringstream str;
        str << line;
        str >> weight_tmp;
        weight_array.push_back(weight_tmp);
    }

    weight = ublas::identity_matrix<double>(weight_array.size());
    double weightall(0.0);
    for (int i = 0; i < weight_array.size(); ++i)
        weightall += weight_array[i];

    for (int i = 0; i < weight_array.size(); ++i)
        weight(i,i) = weight_array[i] / weightall * weight_array.size();

}

void Input::read_eci(const char* file_name){

    Parse_input parse(file_name);
    parse.assignNeed("cluster_index", cluster_index);
    parse.assignNeed("eci", eci);

}

void Input::read_cluster_function_out(){

    Parse_input parse("cluster_function.out");
    parse.assignNeed("n_type", n_type);

    std::ifstream input("cluster_function.out");
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line) && !input.eof()){
        std::stringstream ss;
        ss << line;
        int serial_cluster_number, cluster_number;
        dmatrix function;
        ss >> serial_cluster_number >> cluster_number >> function;
        cluster_function_array.push_back(function);
        cluster_number_array.push_back(cluster_number);
    }

}
void Input::get_temp_array(Parse_input& input, const int& i_temp){

    if (i_temp == 0){
        double temperature;
        input.assignNeed("temp", temperature);
        temp_array.push_back(temperature);
    }
    else if (i_temp == 1){
        double temp_init,temp_final,temp_mul;
        input.assignNeed("temp_init", temp_init);
        input.assignNeed("temp_final", temp_final);
        input.assign("temp_mul", temp_mul, 0.9);
        if (temp_init < temp_final){
            std::cerr << " temp_init is lower than temp_final." << std::endl;
            exit(8);
        }
        if (temp_mul > 1.0 or temp_mul < 0.0){
            std::cerr << " temp_mul must be a value from 0 to 1 !!" 
                << std::endl;
            exit(8);
        }
        temp_array.push_back(temp_init);
        double temp(temp_init);
        while (temp > temp_final){
            temp *= temp_mul;
            temp_array.push_back(temp);
        }
        temp_array.push_back(temp_final);
    }
    else {
        std::cerr << " i_temp must be 0 or 1." << std::endl;
    }

    std::cout << "  temperature = "; 
    for (int i = 0; i < temp_array.size(); ++i){
        std::cout << temp_array[i] << " ";
    }
    std::cout << std::endl;
}

void Input::read_ewald_input
(Parse_input& input, const int& n_atom, const double& volume){

    double weight_factor, accuracy;
    input.assign("weight_factor", weight_factor, 0.001);
    input.assign("accuracy", accuracy, 1e-12);

    double eta_init = pow ((n_atom * weight_factor
                * pow(pi,3) / volume), 0.333333333333);
    input.assign("eta", eta, eta_init);

    double rmax_init = sqrt(-log(accuracy) / eta );
    double gmax_init =  2.0 * sqrt(eta) * sqrt(-log(accuracy));
    input.assign("rmax", rmax, rmax_init);
    input.assign("gmax", gmax, gmax_init);

}

std::string Input::read_file_foption(int argc, char *argv[]){

    std::string output;

    file_name = NULL;
    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"f:"))!=-1){
        switch( option ){
            case 'f' : file_name = optarg;
                       break;
        }
    }
    if (file_name == NULL){
        std::cerr << " File name of structure is needed. (option -f)" 
            << std::endl;
        exit(8);
    }
    output = file_name;

    return output;
}

void Input::set_weight(Parse_input& parse){

    int iweight;
    parse.assign("iweight", iweight, 0);
    if (iweight == 1){
        read_weight("weight");
    }
    else if (iweight == 0){
        weight = ublas::identity_matrix<double>(energy.size());
    }
    else {
        std::cerr << " error : iweight = 0 or 1 " << std::endl;
        exit(8);
    }
}

const int& Input::get_n_sublattice() const{
    return n_sublattice;
}
const std::vector<int>& Input::get_n_type() const{
    return n_type;
}
const int& Input::get_n_lattice() const{
    return n_lattice; 
}
const double& Input::get_symprec() const{
    return symprec; 
}
const std::vector<int>& Input::get_n_atoms_derivative() const{
    return n_atoms_derivative;
}

const std::vector<imatrix>& Input::get_hnf_array() const{
    return hnf_array; 
}
const std::vector<imatrix>& Input::get_snf_array() const{
    return snf_array; 
}
const std::vector<imatrix>& Input::get_left_array() const{
    return left_array; 
}
const vector2_i& Input::get_label_array() const{
    return label_array; 
}
const int& Input::get_max_cluster() const{
    return max_cluster; 
}
const std::vector<double>& Input::get_trunc_distance() const{
    return trunc_distance; 
}
const dmatrix& Input::get_axis_change() const{
    return axis_change;
}
const vector2_pair& Input::get_unique_cluster_array() const{
    return unique_cluster_array;
}
const std::vector<int>& Input::get_spin_array() const{
    return spin_array;
}
const std::vector<int>& Input::get_spin_structure() const{
    return spin_structure;
}
const char* Input::get_file_name() const{
    return file_name;
}
const dvector& Input::get_energy() const{
    return energy;
}
const dmatrix& Input::get_correlation() const{
    return correlation;
}
const dmatrix& Input::get_weight() const{
    return weight;
}
const std::vector<std::string>& Input::get_energy_structure_index() const{
    return energy_structure_index;
}
const std::vector<int>& Input::get_cluster_index() const{
    return cluster_index;
}
const int& Input::get_nloop() const{
    return nloop;
}
const int& Input::get_n_cluster() const{
    return n_cluster;
}
const std::vector<double>& Input::get_temp_array() const{
    return temp_array;
}
const std::vector<int>& Input::get_base_array() const{
    return base_array;
}
const std::vector<int>& Input::get_cluster_candidate() const{
    return cluster_candidate;
}
const int& Input::get_n_pop() const{
    return n_pop;
}
const int& Input::get_max_iteration() const{
    return max_iteration;
}
const int& Input::get_n_elite() const{
    return n_elite;
}
const double& Input::get_p_mating() const{
    return p_mating;
}
const double& Input::get_p_mutation() const{
    return p_mutation;
}
const dvector& Input::get_eci() const{
    return eci;
}
const int& Input::get_composition_cluster() const{
    return composition_cluster;
}
const std::vector<double>& Input::get_composition_definition() const{
    return composition_definition;
}
const std::vector<double>& Input::get_comp_array() const{
    return comp_array;
}
const int& Input::get_error_check() const{
    return error_check;
}
const std::vector<dmatrix>& Input::get_cluster_function_array() const{
    return cluster_function_array;
}
const std::vector<int>& Input::get_cluster_number_array() const{
    return cluster_number_array;
}
const int& Input::get_n_step_ann() const{
    return n_step_ann;
}
const int& Input::get_n_step_eqv() const{
    return n_step_eqv;
}
const std::string& Input::get_simulation_type() const{
    return simulation_type;
}
const std::string& Input::get_init_structure() const{
    return init_structure;
}
const int& Input::get_i_ewald() const{
    return i_ewald;
}
const double& Input::get_epsilon() const{
    return epsilon;
}
const std::vector<double>& Input::get_mu() const{
    return mu;
}
const std::vector<double>& Input::get_charge() const{
    return charge;
}
const std::vector<double>& Input::get_charge_array() const{
    return charge_array;
}
const vector2_i& Input::get_n_change_array() const{
    return n_change_array;
}
const double& Input::get_eta() const{
    return eta;
}
const double& Input::get_rmax() const{
    return rmax;
}
const double& Input::get_gmax() const{
    return gmax;
}
const std::vector<std::string>& Input::get_structure_index_array() const{
    return structure_index_array;
}
const std::vector<int>& Input::get_group_index_array() const{
    return group_index_array;
}
const double& Input::get_init_potential() const{
    return init_potential;
}
const vector2_d& Input::get_ti_mu_array() const{
    return ti_mu_array;
}
const vector2_d& Input::get_ti_temp_array() const{
    return ti_temp_array;
}
const vector2_d& Input::get_ti_comp_array() const{
    return ti_comp_array;
}
const vector2_d& Input::get_ti_energy_array() const{
    return ti_energy_array;
}

const dmatrix& Input::get_axis() const{
    return axis;
}
const std::vector<int>& Input::get_n_atoms() const{
    return n_atoms;
}
const std::vector<dvector>& Input::get_position() const{
    return position;
}
const std::vector<int>& Input::get_type() const{
    return type;
}
const dmatrix& Input::get_axis_primitive() const{
    return axis_primitive;
}
const std::vector<int>& Input::get_n_atoms_primitive() const{
    return n_atoms_primitive;
}
const std::vector<dvector>& Input::get_position_primitive() const{
    return position_primitive;
}
const std::vector<int>& Input::get_type_primitive() const{
    return type_primitive;
}
const std::vector<dvector>& Input::get_position_simulation() const{
    return position_simulation;
}
const double& Input::get_volume() const{
    return volume;
}
const dmatrix& Input::get_reciprocal_axis() const{
    return reciprocal_axis;
}
const std::string& Input::get_criterion() const{
    return criterion;
}
