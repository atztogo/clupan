/***************************************************************************

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
	
    	Header file for mc_binary.cpp

***************************************************************************/

#ifndef __MC_BINARY
#define __MC_BINARY

#include <iostream>
#include <vector>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "gsl/gsl_rng.h"

#include "input.h"
#include "mc_initialize.h"
#include "mc_common.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;

typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<std::vector<int> > > vector3_i;
typedef std::vector<std::vector<std::vector<std::vector<int> > > > vector4_i;

class MC_binary : public MC_common{

    dvector cmc_diff_correlation
        (std::vector<int>& spin, const std::vector<int>& atoms, 
         const MC_initialize& mc_init);

    public:

    MC_binary(const Input& ip, const MC_initialize& mc_init);
    ~MC_binary();

    void cmc(const MC_initialize& mc_init);
    void cmc_ewald(const Input& ip, const MC_initialize& mc_init);
    void find_sqs(const Input& ip, const MC_initialize& mc_init);

};

#endif
