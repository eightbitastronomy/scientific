/**********************************************************************************
 *  kfinitecore.c - library for Kondo calculations at T>0
 *
 *
 *  Author (pseudonomously): eightbitastronomy (eightbitastronomy@protonmail.com)
 *  Copyrighted by eightbitastronomy, 2018.
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
 ***********************************************************************************/

#include "kfinitecore.h"



extern double quad(double (*f)(double,void *), void *, double, double, double, double, char **);






#ifdef USECASEJ32


double                       /* Because V=0, much of this calculation is useless. */
sce ( double x,              /* It only should have J in the case that JV>0, so that with a return value of J, I can know J&V separately. */
      void * params )        /* So, I'll write sceV0 for the case that JV=0. */
{
  struct Parameters * p = params;
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*(15./8.)*(FERMI(XI(x,p->g)+DELTA(x,p->g,p->J,0.),p->mu,p->T)-FERMI(XI(x,p->g)-DELTA(x,p->g,p->J,0.),p->mu,p->T))/DELTA(x,p->g,p->J,0.);
  double gpx = .5*(p->g+x);
  double gmx = .5*(p->g-x);
  double rad = sqrt( gmx*gmx + 14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  return dirac(x)*1.875*( FERMI(gpx-rad,p->mu,KB,p->T) - FERMI(gpx+rad,p->mu,KB,p->T) ) / rad;
}

double                         /* Because V=0, much of this calculation is useless. */
sce_k ( double x,              /* It only should have J in the case that JV>0, so that with a return value of J, I can know J&V separately. */
	double y,              /* So, I'll write sceV0 for the case that JV=0. */
	double z,
	void * params )        
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*(15./8.)*(FERMI(XI(x,p->g)+DELTA(x,p->g,p->J,0.),p->mu,p->T)-FERMI(XI(x,p->g)-DELTA(x,p->g,p->J,0.),p->mu,p->T))/DELTA(x,p->g,p->J,0.);
  double gpn = .5*(p->g+nk);
  double gmn = .5*(p->g-nk);
  double rad = sqrt( gmn*gmn + 14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  return 1.875*( FERMI(gpn-rad,p->mu,KB,p->T) - FERMI(gpn+rad,p->mu,KB,p->T) ) / rad;
}


double
sceV0 ( double x,             /* for the V=0 case */
	void * params )
{
  struct Parameters * p = params;
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*3.75*(FERMI((p->g)/2.,p->mu,p->T)-FERMI(x/2.,p->mu,p->T))/(p->g-x);
  return dirac(x)*3.75*( FERMI(x,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) / (p->g-x);
}

double
sceV0_k ( double x,             /* for the V=0 case */
	  double y,
	  double z,
	  void * params )
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*3.75*(FERMI((p->g)/2.,p->mu,p->T)-FERMI(x/2.,p->mu,p->T))/(p->g-x);
  return 3.75*( FERMI(nk,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) / (p->g-nk);
}

double
sceV0_k_eta( double x,
	     double y,
	     double z,
	     void * params )
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
  double gmn = p->g-nk;
#  ifdef DEBUG
  if (p->flag == 1)
  {
    fprintf(p->fout,"Processing: %f, %f::: %1.8e\n",p->mu,p->T,gmn);
    p->flag = 0;
  }
#  endif
  return 3.75 * ( FERMI(nk,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) * gmn / ( gmn*gmn + (p->eta)*(p->eta) );
}



#else /*S12*/


double
sce( double x,             /* see comment for sce() to see why scespinhalfV0 exists */
     void * params )
{
  struct Parameters * p = params;
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*(.3750)*(FERMI(XI(x,p->g)+DELTAH(x,p->g,p->J,0.),p->mu,p->T)-FERMI(XI(x,p->g)-DELTAH(x,p->g,p->J,0.),p->mu,p->T))/DELTAH(x,p->g,p->J,0.);
  double gpx = .5*(p->g+x);
  double gmx = .5*(p->g-x);
  double rad = sqrt( gmx*gmx + .5625*(p->J)*(p->J)*(p->V)*(p->V) );
  return dirac(x)*.375*( FERMI(gpx-rad,p->mu,KB,p->T) - FERMI(gpx+rad,p->mu,KB,p->T) ) / rad;
}

double
sce_k( double x,             /* see comment for sce() to see why scespinhalfV0 exists */
       double y,
       double z,
       void * params )
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  //return dirac(x)*(.3750)*(FERMI(XI(x,p->g)+DELTAH(x,p->g,p->J,0.),p->mu,p->T)-FERMI(XI(x,p->g)-DELTAH(x,p->g,p->J,0.),p->mu,p->T))/DELTAH(x,p->g,p->J,0.);
  double gpn = .5*(p->g+nk);
  double gmn = .5*(p->g-nk);
  double rad = sqrt( gmn*gmn + .5625*(p->J)*(p->J)*(p->V)*(p->V) );
  return .375*( FERMI(gpn-rad,p->mu,KB,p->T) - FERMI(gpn+rad,p->mu,KB,p->T) ) / rad;
}



double
sceV0( double x,
       void * params )
{
  struct Parameters * p = params;
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  /*return dirac(x)*.750*(FERMI((p->g),p->mu,p->T)-FERMI(x,p->mu,p->T))/(p->g-x);*/
  return dirac(x)*.75*( FERMI(x,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) / (p->g-x);
}

double
sceV0_k( double x,
	 double y,
	 double z,
	 void * params )
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
#  ifdef DEBUG
  if (p->flag == 1)
    {
      fprintf(p->fout,"Processing: %f, %f, %1.10f\n",p->mu,p->T,p->nc);
      p->flag = 0;
    }
#  endif
  return .75*( FERMI(nk,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) / (p->g-nk);
}

double
sceV0_k_eta( double x,
	     double y,
	     double z,
	     void * params )
{
  struct Parameters * p = params;
  double nk = EN(x,y,z);
  double gmn = p->g-nk;
#  ifdef DEBUG
  if (p->flag == 1)
  {
    fprintf(p->fout,"Processing: %f, %f::: %1.8e\n",p->mu,p->T,gmn);
    p->flag = 0;
  }
#  endif
  return .75 * ( FERMI(nk,p->mu,KB,p->T) - FERMI(p->g,p->mu,KB,p->T) ) * gmn / ( gmn*gmn + (p->eta)*(p->eta) );
}


#endif /*USECASEJ32*/





#ifdef USECASEJ32



/* I have not checked this function in a very long time. It is most likely flawed. */

/* ALTERATION: the Jj(j+1)V^2 term should be subtracted (for J<0 condition) not added. DONE, 17FEB2017 */ 
/* SEE BELOW ALSO....all changes have been made, original code left in comments 17feb2017 */
/* NONE OF THESE CHANGES MADE ANY DIFFERENCE. 'DIFF' PERFORMED ON OUTPUTS OF CORRECTED KPMSRHSMOD VS OLD KPMSRHSMOD RETURNED NOTHING */

#define XI  (0.)
#define DELTA(x,y,z,w)  (0.)

double
free_energy( double x,
	     void * params )
{
  struct Parameters * p = params;
  double energy, occfraction, buf_S[4], buf_EN[2];
  /* lambda-plus terms: */
  energy = XI(x,p->g)+DELTA(x,p->g,p->J,p->V);
  occfraction = FERMI(energy,p->mu,KB,p->T);
  if ( occfraction==0. )
    {
      /* it doesn't make sense to perform this calculation... energy*0==0...
	 buf_EN[0]=energy*occfraction; */
      buf_EN[0]=0.;
      buf_S[0]=0.;
    }
  else
    {
      buf_EN[0]=energy*occfraction;
      buf_S[0]=occfraction*log(occfraction);
    }  
  if ( (1.-occfraction)==0. )
    {
      buf_S[1]=0.;
    }
  else
    {
      /* SHOULDN'T THIS READ (1.-occfraction)*log(1.-occfraction) ??? */
      /* buf_S[1]=occfraction*log(occfraction); */
      buf_S[1]=(1.-occfraction)*log(1.-occfraction);
    }
  /* lambda-minus terms: */
  energy = XI(x,p->g)-DELTA(x,p->g,p->J,p->V);
  occfraction = FERMI(energy,p->mu,KB,p->T);
  if ( occfraction==0. )
    {
      /* again, this doesn't make sense. It's 0,
	 buf_EN[2]=energy*occfraction; */
      buf_EN[2]=0.;
      buf_S[2]=0.;
    }
  else
    {
      buf_EN[2]=energy*occfraction;
      buf_S[2]=occfraction*log(occfraction);
    }  
  if ( (1.-occfraction)==0. )
    {
      buf_S[3]=0.;
    }
  else
    {
      /* SHOULDN'T THIS READ (1.-occfraction)*log(1.-occfraction) ???
	 buf_S[3]=occfraction*log(occfraction); */
      buf_S[3]=(1.-occfraction)*log(1.-occfraction);
    }  

  /* It makes a difference, below, whether I have placed JV into J and set V=1, or used a properly calculated J & V from JV,
     entirely because V is squared.
     However, I so far am only using cases where V=0. In these cases, it doesn't make a difference. The E_0 term=0.*/
#ifdef AFFINE
  //return mesh(x)*( sumd(buf_EN,2) - p->J*p->V*p->V + KB*p->T*sumd(buf_S,4) );
  return dirac(x)*( sumd(buf_EN,2) - p->J*p->V*p->V + KB*p->T*sumd(buf_S,4) );
#else
  //return mesh(x)*( sumd(buf_EN,2) - 3.75*p->J*p->V*p->V + KB*p->T*sumd(buf_S,4) );
  return dirac(x)*( sumd(buf_EN,2) - 3.75*p->J*p->V*p->V + KB*p->T*sumd(buf_S,4) );
#endif
}


#else /*S12*/


double
free_energy( double x,
	     void * params )
{
  return 0.;
}


#endif /*USECASEJ32*/




double
fermi( const double xx,
       const double mm,
       const double TT   )
{
  return (1./(exp((xx-mm)/(KB*TT))+1.));
}




double
sumd( double * a, int n )
{
  int sz = sizeof(double);
  double buf = *a;
  while (n>1)
    {
      a+=sz;
      buf+=*a;
      n--;
    }
}
