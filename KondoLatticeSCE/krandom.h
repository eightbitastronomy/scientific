/********************************************************************************************************************************
 *
 *  krandom library: library to provide pseudorandom functions; to hide access to mersenne twister from dSFMT library.
 *
 *  Author: Miguel Abele
 *  KondoLatticeSCE is copyrighted by Miguel Abele, 2018.
 *  dSFMT code is copyrighted by its owners.
 * 
 *  License information:
 *
 *  dSFMT is licensed separately from KondoLatticeSCE. See dSFMT/LICENSE.txt.
 *
 *  This file is a part of KondoLatticeSCE.
 *
 *  KondoLatticeSCE is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 3
 *  of the License, or (at your option) any later version.
 *
 *  KondoLatticeSCE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 ***********************************************************************************
 *
 *  2018 aug 23: moved "shotgun" functions from aug23 version of kondo_bulk to here. Added "multiplier" argument to _slug.
 *
 ********************************************************************************************************************************/


#ifndef KRANDOM_H
#define KRANDOM_H



#ifndef DSFMT_MEXP
#  define DSFMT_MEXP 521
#endif


void start_krand();


void generate_q( double * x,
		 double * y,
		 double * z );

int generate_q_partial( double * x,
			double * y,
			double * z );

double generate_unit();

double shotgun_shot( double a,
		     double b );

double shotgun_slug( double a,
		     double b,
		     double multiplier );


/* unused */
void stop_krand();



#endif /*KRANDOM_H*/
