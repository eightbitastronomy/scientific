/***********************************************************************************************************************
 *
 *  kdefs.h - General definitions.
 *  
 *  Author: Miguel Abele
 *  Copyrighted by Miguel Abele, 2018.
 * 
 *  License information: 
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
 ***********************************************************************************************************************
 *  2018 aug 16: added DEBUG3 into DEBUG decision tree.
 *  2018 jul 27: Placed decision tree: if DEBUG1 or DEBUG2 is defined, DEBUG is then defined so that "fout" member
 *               will also be defined in struct Parameters.
 *  2018 jul 24: All code must now specifiy "USECASEJ32" or not for S=3/2 or S=1/2 respectively, regardless
 *               of whether "ANISOTROPIC" is defined. Previously was only pertinent to ANISOTROPIC, but because
 *               I forgot this, I ran MANY kz's and sz's for S=1/2 using incorrect M and NIMP.
 *  2018 july 4: Extended [LB,UB] energy bounds to [-7.,7], after kz run with J=2.0 demonstrated a need.
 *
 **********************************************************************************************************************/



#ifndef KDEFS

#define KDEFS


#include <stdlib.h>
#include <stdio.h>

#define EN(xx,yy,zz)              (-2.*(cos(xx)+cos(yy)+cos(zz)))
#define ENH(hh,kx,ky,kz)          (-2.*hh*(cos(kx)+cos(ky)+cos(kz)))
#define FERMI(xx,mm,kk,TT)        (1./(exp((xx-mm)/(kk*TT))+1.))
#define FERMIQ(xx,mm,kk,TT)       (1./(expq((xx-mm)/(kk*TT))+1.))
#define HEAVISIDE(xx,yy)          ( ((xx-yy)>=0.) ? 1. : 0. )
#define DOS_WRAPPER(xx)           mesh(xx)
#define NEAR(xx,yy,pp)            ( fabs(xx-yy)<pp )
//#define NEARZERO(xx)              ( (fabs(xx)<ZEROPRECISION) ? (int)1 : (int)0 )
#define NEARZERO(xx)              ( fabs(xx)<ZEROPRECISION )
#define BETWEEN(xx)               ( (xx<1.0-ZEROPRECISION) && (xx>ZEROPRECISION) )
#define SIGNCHANGE(aa,bb)         ( ((signbit(aa)-signbit(bb)) == 0) ? (int)0 : (int)1 )

/* control */
#define REL_ERROR                 (1.e-06)     /* maximum fractional difference, used in numerical integration */
#define PREC                      (1.e-06)     /* precision */
#define NINTINT                   (1000)       /* number of intervals for numerical integration */
#define SUMPREC                   ( .025)       /* interval size for k-summations */
#define MAXSZ                     (256)        /* for re-sizing */
#define ZEROPRECISION             (1e-14)      /* value==zero detection precision */
#define SUMMOR                    k_simpsons_c  /* k-summation function */

/* eV values */
#define KB                        (8.6173e-05)  /* Boltzmann constant */
#define HOPPING                   (1.0000)      /* 't' value, or, hopping coefficient */

/* Other useful parameters */
#define LOWERINT                  (-50.)    /* limit for integration */
#define LOWER                     (-22.)    /* lower bound for rho_m, rho_c integration */
#define UPPER                     (22.)     /* upper bound for rho_m, rho_c integration */
#define LOWERSCE                  (-6.)     /* lower bound for SCE integration */
#define UPPERSCE                  (6.)      /* upper bound for SCE integration */
#define LB                        (-9.0)    /* lower bound for bare conduction integration - INCREASE BOUNDS IF YOU NEED nc's NEARER EMPTY or FULL */
#define UB                        (9.0)     /* upper bound for bare conduction integration - INCREASE BOUNDS IF YOU NEED nc's NEARER EMPTY or FULL */
#define JLOWER                    (-1000.)  /* lowest J value to be tried */

/* S=1/2 vs J=3/2 */
#ifdef ANISOTROPIC
#  ifdef USECASEJ32
#    define NIMP                    ( .25 )        /* impurity occupation */
#    ifdef USEMAGCASE32
#      define MAGM                  (0.250)        /* anisotropic impurity magnetization */
#    else
#      define MAGM                  (0.250)
#    endif /*USEMAGCASE32*/
#    define XI(xx,gg)               (.500)*(gg+xx)
#    define XI_P(xx,gg)             (.500)*(gg-xx)
#    ifdef AFFINE
#      define DELTA(xx,gg,jj,VV)    sqrt(XI_P(xx,gg)*XI_P(xx,gg)+((jj*jj)*(VV*VV)))
#    else
#      define DELTA(xx,gg,jj,VV)    sqrt(XI_P(xx,gg)*XI_P(xx,gg)+((jj*jj)*(14.0625)*(VV*VV)))
#    endif /*AFFINE*/
#  else /* S12 */
#    define NIMP                    ( .5 )
#    define MAGM                    (0.500)
#  endif
#else /*ISOTROPIC*/
#  ifdef USECASEJ32
#    define NIMP                      (.25)
#    define MAGM                      (1.5)
#  else
#    define NIMP                      (0.50)
#    define MAGM                      (0.50)
#  endif
#endif /*ANISOTROPIC*/

/* Debugging decision tree (shrub, really) */
#ifdef DEBUG1
#  define DEBUG
#endif
#ifdef DEBUG2
#  define DEBUG
#endif
#ifdef DEBUG3
#  define DEBUG
#endif

struct Parameters {
  /*mu and J must be first and in this order!!! ...WRT sm.h library and older versions of dispatch... */
  double mu;    /* chemical potential */
  double J;     /* coupling constant */
  double g;     /* impurity kinetic energy */
  double T;
  double nc;
  double V;     /* hybridization */
  double magS;     /* mag is S or M, depending on cond or impurity case, resp. */
  double magM;
  double qx;
  double qy;
  double qz;
  double eta;
#ifdef DEBUG
  FILE * fout;
  int flag;
#endif
};


struct Pwrapper {
  double e;
  struct Parameters * p;
};


/* msolve searching */
enum { LOW,
       MID,
       HI    };





#endif /* KDEFS */
