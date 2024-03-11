/**********************************************************************************
 *  kondo_finite_general.c - Kondo phase calculations at T>0
 *
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
 *******************************************************************************************************************************
 * 22 jun 2018, added tfill_verylow
 *
 * Altered 2018 june 13 as per PSR and sign of linearization constants, etc.
 *  zkpms2.c for zkpms2: revision of zkpms.
 *  Mu-solvers have been changed out for msplit and msolve2, which better handle the flat/plateau regions of mu solutions.
 *  Other small changes are present, such as use of a separate lower integration limit from the lower search limit on mu.
 *  Use of lower and upper mu values from msplit is limited, so far, to choosing the lower mu always.
 *  3 jan 2017
 *  Updates: 10 jan 2017, cleanup of sceVx
 *           11 jan 2017, debug file is created only if DEBUG flag is defined
 *******************************************************************************************************************************/


#include <stdio.h>
#include <stdarg.h>
#include "nint.h"
#include "fio.h"
#include "sm.h"
#include "kfinitecore.h"

/* unistd.h is only used in this file if debugging is needed */
#ifdef DEBUG
#  include <sys/unistd.h>
#endif


#define LOWERMU       ( -5.9999 )
#define INCMU         ( .2 )
#define COUNTMU       ( 61 )
#define DATASZ        ( 1 ) //4
#define RESULTSZ      ( 2 )
#define TFILLER       Tfill_verylow


const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";

enum {
  N,
  R
};

void * masterSCE( void *, void *, void *);
double eta_finder( void * par, double inc );
int Tfill_verylow( double ** d );
int Tfill_low( double ** );
int Tfill_high( double ** );
void mufill(  double ** d, int sz, double lo, double hi, double inc );
void dummyhook(int i);
void egress(int, ... );




int
main (void)
{

  /* Declarations */
  
  int i = 0,j=0,k,l;
  int ret,retval, count;
  char * errorcode = NULL;                     /* for GSL error-codes passed through "quad" function */

  double lowerT, incT ;

  double * tpts;                               /* input data - temperatures */
  double * mupts;                              /* input data - mu / chemical potentials */
  double * results;                            /* return results array */
  double * meshx, * meshy;                     /* dos mesh containers */
  int meshn;                                   /* mesh size */
    
  struct Parameters p;                    /* container for model parameters & values, ktempdefs.h */
  struct d_Function df;                        /* dispatchv2 struct, from sm library */

#ifdef DEBUG
  p.fout = stdout;
  p.flag = 1;
#endif
  

  /* Informational printouts */

#ifdef USECASEJ32
  printf("#J=3/2\n");
  fprintf(stderr,"J=3/2\n");
#  ifdef USEMAGCASE32
  printf("#(Mag 3/2 is defined but useless here)\n");
  fprintf(stderr,"(Mag 3/2 is defined but useless here)\n");
#  endif
#else
  printf("#S=1/2\n");
  fprintf(stderr,"S=1/2\n");
#endif
#ifdef KSUM
  printf("#Summation over k\n");
  fprintf(stderr,"Summation over k\n");
#else
  printf("#Integration\n");
  fprintf(stderr,"Integration\n");
#endif
  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);

  printf("#KPMS V=0 T>0\n#----------------------------------------------------------------------------------------------\n## T        mu           nc          J    \n");


  /* fetch the d.o.s. interpolation mesh */

  fprintf(stderr,"Loading interpolation mesh... ");
  meshn = get_mesh(dirname,meshname,&tpts);
  if ( meshn<1 )
    {
      printf("Mesh problem: %d\n",meshn);
      if (tpts)
	free(tpts);
    }
  ret = mesh_copy( &meshx, &meshy, tpts, meshn );
  if ( ret < 0 )
    {
      if (tpts)
	free(tpts);
    }
  if (tpts)
    {
      free(tpts);
      tpts = NULL;
    }


  /* prepare numerical integration */

  fprintf(stderr,"Preparing numerical integration... ");
  start_ninterp(meshn,meshx,meshy);
  meshx=NULL;
  meshy=NULL;
  start_nint(NINTINT);
  fprintf(stderr,"Done.\n");


  /* Data Prep */
  
  p.V = 0.;
#ifdef KSUM
  //p.eta = eta_finder(&tp,SUMPREC);
#endif
  
  df.f = masterSCE;
  df.param = &p;


  /* Data prep, input arrays */  
  
  fprintf(stderr,"Preparing input data... ");
  count = TFILLER (&tpts) ;
  tpts = malloc(count*sizeof(double)*DATASZ);
  if ( tpts==NULL )
    {
      printf("Failed malloc for temperature points\n");
      egress(1,-1);
    } 
  TFILLER(&tpts);
  mupts = malloc(COUNTMU*sizeof(double));
  if ( tpts==NULL )
    {
      printf("Failed malloc for mu points\n");
      egress(1,-1);
    } 
  for (i=0;i<COUNTMU;i++)
    {
      mupts[i]=LOWERMU+((double)i)*INCMU;
    }
  fprintf(stderr,"Done.\n");
  
  
  /* Data prep & Dispatch loop: */

  fprintf(stderr,"Dispatching with %d (T) x %d (mu) = %d points...",count,COUNTMU,count*COUNTMU);

  for (i=0;i<count;i++)
    {
      p.T = tpts[i];
      
      fflush(stdout);
      fflush(stderr);
      
      dispatchv2(df,COUNTMU,mupts,sizeof(double),&results,RESULTSZ*sizeof(double),egress);
      
      
      /* print results */
      
      for (j=0;j<COUNTMU;j++)
	{
	  if ( isnan(results[j*RESULTSZ+1]) )
	    printf("%f %f %1.8e 00.00\n",tpts[i],mupts[j],results[j*RESULTSZ]);
	  else
	    printf("%f %f %1.8e %1.8e\n",tpts[i],mupts[j],results[j*RESULTSZ],1./results[j*RESULTSZ+1]);
	}
      fflush(stdout);
      
      
      /* clear return-value array, increment, and go */
      
      if (results)
	{
	  free(results);
	  results=NULL;
	}
      
    }
      
  lowerT = 1.;
  incT   = .1;
  
  
  /* clean up */

  fprintf(stderr,"Done. Freeing memory...\n");
  if (tpts)
    {
      free(tpts);
      tpts=NULL;
    }
  if (mupts)
    {
      free(mupts);
      mupts=NULL;
    }
  stop_nint();
  stop_ninterp();

  return EXIT_SUCCESS;
  
};






#ifdef KSUM



void *
masterSCE(  void * params,
	    void * indata,
	    void * outdata    )
{
  char * eecode = NULL;
  struct Parameters * p = params;
  double * input = (double *)indata;
  double * output = (double *)outdata;
  double mu;
  int i=0;
  p->mu = *input;
#  ifdef USECASEJ32
  p->g = p->mu + KB*(p->T)*1.0986122887; /*+kTlog3*/
#  else
  p->g = p->mu;
#  endif
  output[N] = SUMMOR(occupation_bare_finite_k,p,SUMPREC);
  output[R] = SUMMOR(sceV0_k,p,SUMPREC);
#ifdef DEBUG
  fprintf(p->fout,"master: %f %1.6e %1.8e\n",p->T,p->mu,output[N]);
  p->flag = 1;
#endif
  //output[R] = SUMMOR(sceV0_k_eta,p,SUMPREC);
  return output;
}


double eta_finder( void * par, double inc )
{
  struct Parameters * p = par;
  double min = 2*M_PI;
  /* find smallest inc */
  double kx=-M_PI;
  double ky=-M_PI;
  double kz=-M_PI;
  double upperbound=M_PI; /* if i use M_PI+inc, then the integral has the same value twice due to periodicity */
  double buff;
  int i;
  do
    {
      ky=-M_PI;
      do
	{
	  kz=-M_PI;
	  do
	    {
	      buff = fabs(EN(kx,ky,kz) - EN(kx+p->qx,ky+p->qy,kz+p->qz));
	      if ( buff < min )
		{
		  min = buff;
		}
	      kz+=inc;
	    } while (kz<upperbound);
	  ky+=inc;
	} while (ky<upperbound);
      kx+=inc;
    } while (kx<upperbound);

  return min/2.;
}



#else /*integral*/




void *
masterSCE( void * params,
	   void * indata,
	   void * outdata    )
{
  char * eecode = NULL;
  struct Parameters * p = params;
  double * input = (double *)indata;
  double * output = (double *)outdata;
  double mu;
  int i=0;
  p->mu = *input;
#  ifdef USECASEJ32
  p->g = p->mu + KB*(p->T)*1.0986122887; /*+kTlog3*/
#  else
  p->g = p->mu;
#  endif
  output[N] = quad(occupation_bare_finite,p,LOWERSCE,UPPERSCE,0.,REL_ERROR,&eecode);
  if (eecode)
    {
      free(eecode);
      eecode=NULL;
    }
  output[R] = quad(sceV0,p,LOWERSCE,UPPERSCE,0.,REL_ERROR,&eecode);
  if (eecode)
    {
      free(eecode);
      eecode=NULL;
    }
  return output;
}



#endif /*KSUM*/





int Tfill_verylow( double ** d )
{
  int i=0,j=0;
  double T=.001;
  double incT=.001;
  int count = 1000;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  return j;
}




int Tfill_low( double ** d )
{
  int i=0,j=0;
  double T=.01;
  double incT=.01;
  int count = 9;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=.1;
  incT=.02;
  count=45;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=1.00;
  incT=.1;
  count = 90;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=10.0;
  incT=.5;
  count = 180;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=100.;
  incT=2.;
  count = 450;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  return j;
}


int Tfill_high( double ** d )
{
  int i=0,j=0;
  double T=.02;
  double incT=.02;
  int count = 5;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=.1;
  incT=.1;
  count=10;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=1.00;
  incT=.5;
  count = 18;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=10.;
  incT=1.;
  count = 90;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=100.;
  incT=5.;
  count = 180;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=1000.;
  incT=50.;
  count = 180;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  T=10000.;
  incT=500.;
  count = 180;
  for ( i=0; i<count; i++ )
    {
      if ((*d)!=NULL)
	(*d)[j*DATASZ]=T+((double)i)*incT;
      j++;
    }
  return j;
}

void mufill(  double ** d, int sz, double lo, double hi, double inc )
{
  int j;
  if ((*d)!=NULL)
    {
      for (j=0;j<sz;j++)
	{
	  (*d)[j*DATASZ+1]=lo;
	  (*d)[j*DATASZ+2]=hi;
	  (*d)[j*DATASZ+3]=inc;
	}
    }
}


void egress(int ac, ... )
{
  va_list alist;
  int i;
  int retval;
  void * container;
  if ( ac==0 )
    return;
  va_start(alist,ac);
  /* first list-argument must be the return value from dispatch */
  retval = va_arg(alist,int);
  /* all subsequent arguments are pointers to memory which must be freed */
  for (i=1;i<ac;i++)
    {
      container = (void *)va_arg(alist,void *);
      if ( container )
	free(container);
    }
  if ( retval == 0 )
    return;
  else
    {
      finish_files();
      stop_nint();
      stop_ninterp();
      exit(retval);
    }
}
