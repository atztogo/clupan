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

	    Class for selecting structures
	
*****************************************************************************/

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include "input.h"

namespace ublas = boost::numeric::ublas;

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc, char *argv[]){

    const char *file_name = NULL;
    int n_select(0), n_index(0);

    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"g:r:n:"))!=-1){
        switch( option ){
            case 'g' : file_name = optarg;
                       break;
            case 'r' : n_index = atoi(optarg);
                       break;
            case 'n' : n_select = atoi(optarg);
                       break;
        }
    }
    if (file_name == NULL and n_index < 1){
        std::cerr << " Specify file name of group (option -g) "
            << "or structure index for random selection (option -r)." 
            << std::endl;
        exit(8);
    }
    if (n_select < 1){
        std::cerr << " n_select should be a positive value (option -n)." 
            << std::endl;
        exit(8);
    }

    std::vector<std::string> all_index;
    if (file_name != NULL){

        Input ip;
        ip.read_group(file_name);

        std::vector<std::string> structure = ip.get_structure_index_array();
        std::vector<int> group = ip.get_group_index_array();

        int max_group = *max_element(group.begin(), group.end());
        std::vector<std::vector<std::string> > structure_group(max_group);
        for (int i = 0; i < group.size(); ++i){
            structure_group[group[i]-1].push_back(structure[i]);
        }

        srand(clock()+getpid());

        for (int i = 0; i < structure_group.size(); ++i){
            std::cout << "  number of structures in group " << i+1 
                << " = " << structure_group[i].size() << std::endl;
        }
        for (int i = 0; i < max_group; ++i){
            std::set<int> rand_index;
            while (rand_index.size() != n_select){
                int index = rand() % int(structure_group[i].size());
                rand_index.insert(index);
            }
            std::set<int>::iterator it = rand_index.begin();
            while (it != rand_index.end()){
                all_index.push_back(structure_group[i][*it]);
                ++it;
            }
        }
    }
    else if (n_index > 0){
        srand(time(NULL));
        std::set<std::string> rand_index;
        while (rand_index.size() != n_select){
            std::stringstream ss;
            std::string output;
            int index = rand() % int(n_index) + 1;
            ss << index;
            ss >> output;
            rand_index.insert(output);
        }
        all_index.resize(n_select);
        std::copy(rand_index.begin(), rand_index.end(), all_index.begin());
    }

    std::cout << "  randomly selected structures = ";
    for (int i = 0; i < all_index.size(); ++i){
        if (i == all_index.size() - 1 ){
            std::cout << all_index[i];
        }
        else {
            std::cout << all_index[i] << ",";
        }
    }
    std::cout << std::endl;

    return 0;

}
