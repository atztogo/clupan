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

	    Header file for ga.cpp

******************************************************************************/

#ifndef __GA
#define __GA

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include <gsl/gsl_rng.h>

#include "least_squares.h"
#include "cv.h"

typedef std::vector<std::vector<int> > vector2_i;

class GA{

	public: 

	GA();
	~GA();

    std::vector<int> initial_combination
        (const std::vector<int>& candidate, const int& n_select, gsl_rng * r);

    double calc_score(const dvector& energy, 
        const dmatrix& correlation_all, const dmatrix& weight, 
        const std::vector<int>& base_array, 
        const std::vector<int>& combination);

    std::vector<int> change_combination
        (const std::vector<int>& candidate, 
         const std::vector<int>& combination_old, gsl_rng * r);

    std::vector<int> generate_child
        (const std::vector<int>& parent1, const std::vector<int>& parent2,
         const std::vector<int>& candidate, const double& p_mating, 
         const double& p_mutation, const int& n_elite, gsl_rng * r);

};

#endif
