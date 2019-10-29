/**********************************************************************************
 *  afm_bulk.c - AFM phase bulk/nonboundary calculations
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
 *********************************************************************************************************************************
 *
 *  2018sep13: Changed name 'CRUDELUMP' to 'LUMP'.
 *  2018sep7:  Added CRUDELUMP. All test-data is added to to-do list, not just one T at a time. Hence, CRUDELUMP does not allow
 *             for basing next guess upon previous result. This will allow better use of nodes on Owl's nest.
 *  2018sep4:  Modifications to bulkcore & kondo_bulk on aug23 left this code broken, bc here mu_finder is expected to modify
 *             p->nc in place. Altered this code to recalculate nc.
 *  2018aug20: dump_failures modified to print a newline character.
 *  2018aug17: added dump_failures. Waiting on a test run.
 *  2018aug16: corrections made to free energy calls & functions.
 *  2018aug13: free energy calculation altered in bulkcore.c and thus altered in solver, here.
 *  2018aug9:  corrected the temperature incrementer near end of main(). Will ignore LOWT if 1.0, will add to LOWT if not.
 *             (having altered jul31 code on owlsnest to increment without LOWT, I ended up with higher T runs that used, e.g.,
 *             3750,250,500, etc, instead of 3750,4000,4250, etc.)
 *  2018jul28: solver loop criterion is input[MM] ~= p->magM, but input[MM] never changes. 
 *             Altered so to use "last magM" ~= p->magM.
 *  2018jul27: altered the 2nd level of debug output, DEBUG2, in solver. Changed DEBUG -> DEBUG1 in accordance with changes
 *             to kdefs.h.
 *
 *******************************************************************************************************************************/




#include <stdio.h>
#include <stdarg.h>
#include "bulkcore.h"
#include "nint.h"
#include "sm.h"
#include "fio.h"



#define JCOUP              ( 1. )
#define QSPLIT             8      /* must be even; splitting increments are to be calculated out of 2pi, but split is on [0,pi] */
#define COUNTT             1
#define INCT               (250.)
#define LOWT               (1.0)
#define STARTnc            (3.1)
#define STARTmagS          (.8250547)//(.929313)
#define STARTmagM          (1.5)
#define STARTmu            (2.644546)//(2.277)
#define STARTglam          (3.060265)//(.883)
#define ETA                (.01)
#define FAILNUMBER         5
#define DATAFIELD          8
#define RESULTFIELD        6

//T        x        y        z        nc2         mu          g           M           S           E
//1.0000 2.3555 2.3555 1.5702 2.900006e+00 2.276537e+00 2.744377e+00 1.500000e+00 9.293091e-01 -3.534887e-04
//1.0000 1.5702 1.5702 1.5702 3.099999e+00 2.644546e+00 3.060265e+00 1.500000e+00 8.250547e-01 -3.121566e-04

#define MU_FUNCTION        mu_finder


enum {
  T,
  MS,
  MM,
  MU,
  G,
  QX,
  QY,
  QZ
};


enum {
  ONC,
  OMU,
  OG,
  OMS,
  OMM,
  OF
};



const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";


void * solver( void * params, void * inputdat, void * outputres );
void dump_failures( int n, double * in, int infields, double * out, int outfields, int checkfield, FILE * f );
void egress(int, ... );




int
main(void)
{

  struct Parameters p;
  struct d_Function df;
  double * results = NULL;
  double * dpts = NULL;
  int i,j,k,l,m,n;
  double incq = ((2*M_PI)/QSPLIT)-PREC;
  int datacount;
  int record;
  FILE * input, * mdata;
  char * errorcode;

#ifdef DEBUG
  p.fout = stdout;
#endif


  /* Count q-splitting for datacount: count up to qsplit-1, then add one more iteration */

  k=0;
  for (l=0;l<(QSPLIT/2);l++)
    {
      for (m=0;m<=l;m++)
	{
	  for (n=0;n<=m;n++)
	    {
	      k++;
	    }
	}
    }
  datacount = k+1;

#ifdef LUMP
  datacount *= COUNTT ;
#endif  

  /* print banner info */
  
  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
#ifdef LUMP
  printf("#Data count x T count (q-split %d): %d\n",QSPLIT,datacount);
  fprintf(stderr,"Data count x T count (q-split %d): %d\n",QSPLIT,datacount);
#else
  printf("#Count per iteration of T (q-split %d): %d\n",QSPLIT,datacount);
  fprintf(stderr,"Count per iteration of T (q-split %d): %d\n",QSPLIT,datacount);
#endif
  fprintf(stdout,"#J=%f, nc=%f\n##T        x        y        z        nc2         mu          g           M           S           E\n",JCOUP,STARTnc);
  

  /* prepare numerical integration */

  start_nint(NINTINT);


  /* prepare inputs */
  
  dpts = malloc( DATAFIELD*datacount*sizeof(double) );
  if ( dpts==NULL )
    {
      fprintf(stderr,"Error allocating AFM memory\n");
      stop_ninterp();
      stop_nint();
      return EXIT_SUCCESS;
    }

  p.T = LOWT;
  p.nc = STARTnc;
  p.J=JCOUP;
  p.eta = ETA;
  df.f = solver;
  df.param = &p;

  k=0;
#ifdef LUMP
  for (j=0;j<COUNTT;j++)
    {
#endif
      for (l=0;l<(QSPLIT/2);l++)
	{
	  for (m=0;m<=l;m++)
	    {
	      for (n=0;n<=m;n++)
		{
#ifdef LUMP
		  if ( NEAR(LOWT,1.,PREC) )
		    {
		      if ( j==0 )
			dpts[k++] = LOWT;
		      else
			dpts[k++] = ((double)j)*INCT;
		    }
		  else
		    {
		      dpts[k++] = LOWT + ((double)j)*INCT;
		    }		  
#else
		  dpts[k++] = LOWT;
#endif
		  dpts[k++] = STARTmagS;
		  dpts[k++] = STARTmagM;
		  dpts[k++] = STARTmu;
		  dpts[k++] = STARTglam;
		  dpts[k++] = ((double)l)*incq-4.*PREC;
		  dpts[k++] = ((double)m)*incq-4.*PREC;
		  dpts[k++] = ((double)n)*incq-4.*PREC;
		}
	    }
	}
      l=QSPLIT/2;
      m=0;
      n=0;
#ifdef LUMP
      if ( NEAR(LOWT,1.,PREC) )
	{
	  if ( j==0 )
	    dpts[k++] = LOWT;
	  else
	    dpts[k++] = ((double)j)*INCT;
	}
      else
	{
	  dpts[k++] = LOWT + ((double)j)*INCT;
	}
#else
      dpts[k++] = LOWT;
#endif
      dpts[k++] = STARTmagS;
      dpts[k++] = STARTmagM;
      dpts[k++] = STARTmu;
      dpts[k++] = STARTglam;
      dpts[k++] = ((double)l)*incq-4.*PREC;
      dpts[k++] = ((double)m)*incq-4.*PREC;
      dpts[k++] = ((double)n)*incq-4.*PREC;
#ifdef LUMP
    }
#endif
  
  
  /* T loop */

#ifndef LUMP
  i=0;
    
  while ( i<COUNTT )
    {
#endif

      /* flush streams */
      
      fflush(stdout);
      fflush(stderr);
      
      /* dispatch */
      
      dispatchv2( df, datacount, dpts, DATAFIELD*sizeof(double), &results, RESULTFIELD*sizeof(double), egress );

#ifdef LUMP

      for ( l=0;l<datacount;l++ )
	{      /* T      x     y     z     nc2   mu    g     M     S     E*/
	  printf("%1.4f %1.4f %1.4f %1.4f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n",
		 dpts[l*DATAFIELD+T],dpts[l*DATAFIELD+QX],dpts[l*DATAFIELD+QY],dpts[l*DATAFIELD+QZ],
		 results[l*RESULTFIELD+ONC],results[l*RESULTFIELD+OMU],results[l*RESULTFIELD+OG],
		 results[l*RESULTFIELD+OMM],results[l*RESULTFIELD+OMS],results[l*RESULTFIELD+OF]);
	}
      fflush(stdout);

#else
      
      /* dump failures */

      dump_failures( datacount,dpts,DATAFIELD,results,RESULTFIELD,ONC,stdout );

      
      /* minimize F */
      
#  ifdef DEBUG1
      fprintf(p.fout,"Minimizing...\n");
      fflush(p.fout);
#  endif

      record = filter_minimum( results, datacount, RESULTFIELD, OF );
      
      /* dump best */

      /*      T     x     y     z     nc2   mu    g     M     S     E    */
      printf("%1.4f %1.4f %1.4f %1.4f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n",
	     dpts[record*DATAFIELD+T],dpts[record*DATAFIELD+QX],dpts[record*DATAFIELD+QY],dpts[record*DATAFIELD+QZ],
	     results[record*RESULTFIELD+ONC],results[record*RESULTFIELD+OMU],results[record*RESULTFIELD+OG],
	     results[record*RESULTFIELD+OMM],results[record*RESULTFIELD+OMS],results[record*RESULTFIELD+OF]   );
#endif
      
      
      /* set params for next iteration */

#ifdef DEBUG1
      fprintf(p.fout,"Resetting input data...\n");
      fflush(p.fout);
#endif

#ifndef LUMP
      i++;

      for (l=0;l<datacount;l++)
	{
	  if ( NEAR(LOWT,1.,PREC) )
	    {
	      dpts[l*DATAFIELD+T] = ((double)i)*INCT;
	    }
	  else
	    {
	      dpts[l*DATAFIELD+T] = LOWT + ((double)i)*INCT;
	    }
	  dpts[l*DATAFIELD+MS] = results[record*DATAFIELD+OMS];
	  dpts[l*DATAFIELD+MM] = results[record*DATAFIELD+OMM];
	  dpts[l*DATAFIELD+MU] = results[record*DATAFIELD+OMU];
	  dpts[l*DATAFIELD+G] = results[record*DATAFIELD+OG];
	}

      /* clear results */
  
      if (results)
	{
	  free(results);
	  results=NULL;
	}

    }
#endif
  

  /* flush streams */
  
  fflush(stdout);
  
  
  /* clean up */
  
  if (dpts)
    {
      free(dpts);
    }
  if (results)
    {
      free(results);
    }

  stop_nint();

  return EXIT_SUCCESS;

};








void * solver( void * params,
	       void * inputdat,
	       void * outputres )
{
  struct Parameters * p = params;
  double * input = inputdat;
  double * output = outputres;
  double runningnc, lastmm;
  char * ecode = NULL;
  int iter = 0;

  p->T = input[T];
  p->magS = input[MS];
  p->magM = input[MM];
  p->mu = input[MU];
  p->g = input[G];
  p->qx = input[QX];
  p->qy = input[QY];
  p->qz = input[QZ];

#ifdef DEBUG2
  fprintf(p->fout,"--------------------------------------------------------\nT %f, S %f, M %f, mu %f, g %f, (%f,%f,%f)\n",p->T,p->magS,p->magM,p->mu,p->g,p->qx,p->qy,p->qz);
  fflush(p->fout);
#endif
  /* guess MagM: first guess is the input */

  while ( iter < FAILNUMBER )
    {
#ifdef DEBUG2
      fprintf(p->fout,"        %d:\n",iter);
      fflush(p->fout);
#endif
      p->mu = MU_FUNCTION( mu_wrapper_afm,p,STARTnc );
      /* aug23 code, nc is not stored in p. Nc is not calculated exactly for the returned mu, so MUST be recalculated */
      runningnc = quad( occupation_afm_c_wrapper,p,LOWER,UPPER,0,REL_ERROR,&ecode );
      if (ecode)
	{
	  free(ecode);
	  ecode=NULL;
	}
      p->magS = quad( sce_afm_s_wrapper,p,LOWER,UPPER,0,REL_ERROR,&ecode );
      if (ecode)
	{
	  free(ecode);
	  ecode=NULL;
	}
      p->g = glam_finder( glam_wrapper_afm,p,1.0 );   /* nf constraint is solved here */
      lastmm = p->magM;
      p->magM = sce_afm_m( p->g,p->mu,p->J,p->magS,p->T );
#ifdef DEBUG2
      fprintf(p->fout,"        mu %1.4e, nc %1.4e S %1.4e, g %1.4e, M %1.4e\n",p->mu,runningnc,p->magS,p->g,p->magM);
      fflush(p->fout);
#endif
      if ( NEAR(lastmm,p->magM,PREC)==1 && NEAR(runningnc,STARTnc,PREC)==1 )
	{
	  break;
	}
      iter++; /* must be last in loop because below I will check the value of iter */
    }

#ifdef DEBUG2
  fprintf(p->fout,"        ...iter=%d...",iter);
  fflush(p->fout);
#endif

  if ( iter >= FAILNUMBER )
    {
      output[ONC] = -1.;
      output[OMU] = 0.;
      output[OG]  = 0.;
      output[OMS] = -1.;
      output[OMM] = -1.;
      output[OF] = 0.;
    }
  else
    {
      output[ONC] = runningnc;
      output[OMU] = p->mu;
      output[OG]  = p->g;
      output[OMS] = p->magS;
      output[OMM] = p->magM;
      output[OF] = free_energy_afm( p );
    }
#ifdef DEBUG2
  fprintf(p->fout," F %1.4e\n",output[OF]);
  fflush(p->fout);
#endif
  return (void *)output;
}









void dump_failures( int n,
		    double * in,
		    int infields,
		    double * out,
		    int outfields,
		    int checkfield,
		    FILE * f       )
{
  int i,j;
  if ( in==NULL )
    return;
  if ( out==NULL )
    return;
  if ( n<1 )
    return;
  if ( infields<1 )
    return;
  if ( outfields<1 )
    return;
  for ( i=0; i<n; i++ )
    {
      if ( out[i*RESULTFIELD+checkfield] < 0. )
	{
	  for ( j=0;j<DATAFIELD;j++ )
	    {
	      fprintf(f,"%1.4e ",in[i*DATAFIELD+j]);
	    }
	  fprintf(f,"| ");
	  for ( j=0;j<RESULTFIELD;j++ )
	    {
	      fprintf(f,"%1.4e ",out[i*RESULTFIELD+j]);
	    }
	  fprintf(f,"\n");
	}
    }
  return;
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
	{
	  free(container);
	  container=NULL;
	}
    }
  if ( retval == 0 )
    return;
  else
    {
      finish_files();
      stop_nint();
/*
      stop_ninterp();
*/
      exit(retval);
    }
}
