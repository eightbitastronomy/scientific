/**********************************************************************************
 *  afmiso_finite.c - SDW phase calculations at T>0
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
 ***********************************************************************************
 *
 *  Nc is [0,4] !
 *
 *  29 jun 2018, lower/upper limits added to mufinder. At AFM-PM boundary, the limits don't change.
 *  2018 jun 7:
 *    o  KSUM references removed; this code has no scalar-integration option. KSUM removed from afmisofinitecore.c.
 *    o  Code converted to process data in one lump. Use -DLUMP to enable this mode.
 *
 ***********************************************************************************************************************/


#include <stdio.h>
#include <stdarg.h>
#include "afmisofinitecore.h"
#include "nint.h"
#include "sm.h"
#include "fio.h"



#define JCOUP              ( 1. )
#define LOWERNC            ( .68 )
#define INCNC              ( .02 )
#define COUNTNC            ( 2 )
#define QSPLIT             ( 4 )
#define MAXT               ( 10000. )
#define R                  2 /*modifies the sdwcore.h enumerations*/
#define RESULTFIELD        3
#ifdef LUMP
#  define DATAFIELD        4
#else
#  define DATAFIELD        1
#endif



const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";



void * bissection_of_T_n( void * par, void * in, void * out );
void egress(int, ... );




int
main(void)
{

  struct Parameters p;
  struct d_Function df;
  double * results = NULL;
  double * dpts = NULL;
  int i,j,k,l,m;
  double qx,qy,qz;
  double incq = ((2*M_PI)/QSPLIT)-PREC;
  int qtotal = ((QSPLIT*QSPLIT*QSPLIT)+(3*QSPLIT*QSPLIT)+(2*QSPLIT))/6;
  int datacount;
  FILE * input, * mdata;

  
  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  printf("#%d divisions of q: %d calculations\n",QSPLIT,countq(QSPLIT));
  fprintf(stderr,"%d divisions of q: %d calculations\n",QSPLIT,countq(QSPLIT));
  fprintf(stdout,"#J=%f\n##qx      qy      qz      nc      T         mu          sce\n",JCOUP);
  fprintf(stderr,"J=%f, (qx,qy,qz,nc,T,mu,sce)\n",JCOUP);

  
  /* prepare numerical integration */

  start_nint(NINTINT);


  /* prepare input values */

  fprintf(stderr,"Preparing input values...\n");
#ifdef LUMP
  dpts = malloc(qtotal*COUNTNC*DATAFIELD*sizeof(double));
  l=0;
  for (m=0;m<COUNTNC;m++)
    {
      for (i=0;i<QSPLIT;i++)
	{
	  for (j=0;j<=i;j++)
	    {
	      for (k=0;k<=j;k++)
		{
		  dpts[l++]=LOWERNC+((double)m)*INCNC;
		  dpts[l++]=((double)i)*incq-M_PI;
		  dpts[l++]=((double)j)*incq-M_PI;
		  dpts[l++]=((double)k)*incq-M_PI;
		}
	    }
	}
    }
  datacount = qtotal*COUNTNC ;
#else
  dpts = malloc(COUNTNC*DATAFIELD*sizeof(double));
  for (i=0; i<COUNTNC; i++)
    {
      dpts[i]=LOWERNC+INCNC*((double)i);
    }
  datacount = COUNTNC;
#endif
  fprintf(stderr,"Entering loop.\n");


  /* prepare inputs */

  p.J=JCOUP;
  df.f = bissection_of_T_n;
  df.param = &p;
  

  /* flush streams */
  
  fflush(stdout);
  fflush(stderr);


  /* q loop */

#ifndef LUMP
  for (i=0;i<QSPLIT;i++)
    {
      for (j=0;j<=i;j++)
	{
	  for (k=0;k<=j;k++)
	    {

	      p.qx=((double)i)*incq-M_PI;
	      p.qy=((double)j)*incq-M_PI;
	      p.qz=((double)k)*incq-M_PI;
#endif

	      /* dispatch */

	      dispatchv2( df, datacount, dpts, DATAFIELD*sizeof(double), &results, RESULTFIELD*sizeof(double), egress );

	      /* dump results */
	      

	      for (l=0;l<datacount;l++)
		{
		  /*      x     y     z     nc    T     mu    R     */
		  printf("%1.4f %1.4f %1.4f %1.3f %1.6e %1.6e %1.6e\n",
#ifdef LUMP
			 dpts[l*DATAFIELD+1],dpts[l*DATAFIELD+2],dpts[l*DATAFIELD+3],dpts[l*DATAFIELD],
#else
			 p.qx,p.qy,p.qz,dpts[l],
#endif
			 results[l*RESULTFIELD+T],results[l*RESULTFIELD+M],results[l*RESULTFIELD+R]  );
		}
	      
	      /* flush streams */
	      
	      fflush(stdout);
	        
	      /* clear results */

	        if (results)
		  {
		    free(results);
		    results=NULL;
		  }

#ifndef LUMP
	    } /* for k */
	}     /* for j */
    }         /* for i */
#endif

  
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






void * bissection_of_T_n( void * par, void * in, void * out )
{
  struct Parameters * p = par;
  double * retval = out;
  double * input = in;
  double lowerT = PREC;
  double upperT = MAXT;
  double resultSCE;
  double mu;
  char * ecode = NULL;
  double lowerdomain;
#ifndef LUMP
  p->nc = *input;  /* qx,qy,qz are set in main for each fork */
#else
  p->nc = input[0];
  p->qx = input[1];
  p->qy = input[2];
  p->qz = input[3];
#endif
  lowerdomain = LB; /* at the AFM-PM boundary, M ~ 0, so for the band lower-limit, -6-JM = -6. */
  
  /* bissection assuming that RHS crosses 1 only once */
  
  do
    {
      p->T = (upperT+lowerT)/2.;
      mu = mufinder(p,p->nc,lowerdomain,-1.*lowerdomain);
      p->mu = mu;
      resultSCE = magRHS(p,lindhard_boundary);
      if ( resultSCE>1. )
	{
	  lowerT = p->T;
	}
      else
	{
	  upperT = p->T;
	}
    } while ( NEAR(upperT,lowerT,PREC)==0 );
  retval[T] = p->T;
  retval[M] = p->mu;
  retval[R] = resultSCE;
  return (void *) retval;
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
      exit(retval);
    }
}
