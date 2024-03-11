/**********************************************************************************
 *  kondo_bulk.c - Kondo phase bulk/nonboundary calculations
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
 ***************************************************************************************************************************
 *
 *  2018 aug 23: stability modifications to rebuild of aug22. Cleaned up.
 *  2018 aug 22: ongoing rebuild of solver & solver_s functions. Stable output achieved in version "2"'s of solvers.
 *  2018 aug 20: replaced kondo energy function in solver & solver_s
 *  2018 aug 14: dispatch loop altered to be a while loop with a post-dispatch increment of i, so that LOWT determines
 *               the progression thereafter (as in afm_bulk.c). Merged with aug13.
 *  2018 aug 13: altered free_energy calc in bulkcore.c and thus altered it here also
 *  2018 aug 9:output altered to include inputnc
 *  2018 aug 8: solver() updated to include changes found in solver_s2()
 *  2018 jul 26: jul25 update was applied only to solve, not solve_s
 *  2018 jul 25: testing wider bounds for randomizing g/lam and V during solve
 *
 ***************************************************************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include "bulkcore.h"
#include "nint.h"
#include "sm.h"
#include "fio.h"
#include "krandom.h"



#define JCOUP              ( 1. )
#define COUNTT             2
#define INCT               (250.)
#define LOWT               (500.0)
#define STARTnc            (3.214)
#define STARTV             (.5963258)//(.67024)//(.61)(.22)
#define STARTmu            (4.763712)//(3.2513)//(-1.8)
#define STARTglam          (4.196838)//(3.5598)//(-1.5)
#define ETA                (.01)
#define FAILNUMBER         500
#define RETRYNUMBER        6
#define DATAFIELD          4
#define RESULTFIELD        6
#define INPUTFIELD         15
#define UPPERLIMITV        ( 10. )



enum {
  NC,
  V,
  MU,
  G
};


enum {
  ONC,
  ONF,
  OMU,
  OG,
  OV,
  OF
};


enum {
  FG,
  FJV,
  FMU,
  FNF,
  FNC,
  FJ,
  FV,
  FE0,
  FEN,
  FEM,
  FC,
  FEI,
  FEE,
  F4N,
  F4E
};




const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";
                 /*      g   jv  mu  nf  nc  J   V   E0  EN  EM  C   EIG E   4nc 4E       */	  
const char format[61] = "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n\0";


void * solver( void * params,
	       void * inputdat,
	       void * outputres );
int solver_s( struct Parameters * p,
	      double * input,
	      double * output );
double bissect_V( void * param,
		  double lb,
		  double ub       );
int glam_looper( double (*gf)(void *, double),
		 double (*mf)(void *, double),
		 void * par,
		 double nimp );
void egress(int, ... );






int
main( int argc, char ** argv )
{

  struct Parameters p;
#ifndef SINGLEPROC
  struct d_Function df;
#endif
  double * results = NULL;
  double * dpts = NULL;
  double dblinput[INPUTFIELD];
  int i,j,k,l,m,n;
  int ret;
  int datacount;
  int record;
  FILE * input, * mdata;
  char * errorcode;

#ifdef DEBUG
  p.fout = stdout;
#endif


  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stdout,"#J=%f\n##T       nc-in    nc-out      nf          mu           V           g            E\n",JCOUP);


  /* open file from command-line */
  
  if ( argc == 2 )
    {
      fprintf(stderr,"Opening %s...\n",argv[1]);
      input = fopen(argv[1],"r");
      if ( input==NULL )
	{
	  fprintf(stderr,"Error opening file.\n");
	  return EXIT_SUCCESS;
	}
    }
  else
    {
      fprintf(stderr,"Couldn't open a file for input\n");
      return EXIT_SUCCESS;
    }


  /* count lines in input */

  datacount = 0;
  while ( feof(input) != EOF )
    {
      ret = getline(&errorcode,(size_t *)&n,input);
      if ( ret > 0 )
	{
	  datacount++;
	  if ( errorcode )
	    {
	      free(errorcode);
	      errorcode=NULL;
	    }
	}
      else
	break;
    }
  rewind(input);

  
  /* allocate space for input */

  dpts = malloc( DATAFIELD*datacount*sizeof(double) );
  if ( dpts==NULL )
    {
      fprintf(stderr,"Error allocating input memory\n");
      stop_nint();
      return EXIT_SUCCESS;
    }


  /* fetch & prepare the kz data */

  i = 0;
  while ( feof(input) != EOF )
    {
      
      /* Read a line and dispatch if EOF has been reached */
      /*                  g   jv  mu  nf  nc  J   V   E0  EN  EM  C   EIG E   4nc 4E       */	     
      ret = fscanf(input,format,
		   &(dblinput[FG]),&(dblinput[FJV]),&(dblinput[FMU]),
		   &(dblinput[FNF]),&(dblinput[FNC]),&(dblinput[FJ]),
		   &(dblinput[FV]),&(dblinput[FE0]),&(dblinput[FEN]),
		   &(dblinput[FEM]),&(dblinput[FC]),&(dblinput[FEI]),
		   &(dblinput[FEE]),&(dblinput[F4N]),&(dblinput[F4E])  );
      if (ret < 1 )
	{
	  fprintf(stderr,"Done.\n");
	  break;
	}

      dpts[i*DATAFIELD+NC] = dblinput[F4N];  /* not FNC because bulk code is not per-spin-channel */
      dpts[i*DATAFIELD+V] = dblinput[FV];
      dpts[i*DATAFIELD+MU] = dblinput[FMU];
      dpts[i*DATAFIELD+G] = dblinput[FG];
      i++;
    }
  fclose(input);
  

#ifdef SINGLEPROC

  /* prepare output data */
  
  results = malloc( RESULTFIELD*sizeof(double) );
  if ( results==NULL )
    {
      printf("Error allocating memory for results\n");
      if (dpts)
	free(dpts);
      stop_nint();
      return EXIT_FAILURE;
    }
#endif
  
  
  /* prepare input parameters */
  
  p.T = LOWT;
  p.eta = ETA;
  p.J = JCOUP;
#ifndef SINGLEPROC
  df.f = solver;
  df.param = &p;
#endif


  /* prepare numerical integration */

  start_nint(NINTINT);
  start_krand();

  
  /* flush streams */
  
  fflush(stdout);
  fflush(stderr);


  /* Temperature loop */

  i=0;

  while ( i<COUNTT )
    {
  
      /* dispatch */

#ifdef SINGLEPROC
      ret = 0;
      int flag = 0;
      for ( j=0;j<datacount;j++ )
	{
	  if (flag)
	    {
	      dpts[j*DATAFIELD+MU] = 1.574238;
	      dpts[j*DATAFIELD+V] = .6290369;
	      dpts[j*DATAFIELD+G] = 3.195885;
	      flag = 0;
	    }
	  if ( (j>0)&&(ret==0) )
	    { 
	      dpts[j*DATAFIELD+MU] = results[OMU];
	      dpts[j*DATAFIELD+V] = results[OV];
	      dpts[j*DATAFIELD+G] = results[OG];
	    }
	  fprintf(stderr,"%1.4f %1.6f ------------ %1.6e %1.6e %1.6e --------- :\n",
		  p.T, dpts[j*DATAFIELD+NC],
		  dpts[j*DATAFIELD+MU],dpts[j*DATAFIELD+V],dpts[j*DATAFIELD+G] );
	  ret = solver_s( &p,&(dpts[j*DATAFIELD]),results );
	  /*      T      nc    nf    mu    V     g     E    */
	  printf("%1.4f %1.6f %1.6e %1.6e %1.6e %1.6e %1.6e\n",
		 p.T,
		 results[ONC],results[ONF],results[OMU],
		 results[OV],results[OG],results[OF]  );
	  fprintf(stderr,"%1.4f %1.6f %1.6e %1.6e %1.6e %1.6e %1.6e\n",
		  p.T,
		  results[ONC],results[ONF],
		  results[OMU],results[OV],results[OG],results[OF]  );
	  fflush(stdout);
	  fflush(stderr);
	}

#else      

      dispatchv2( df, datacount, dpts, DATAFIELD*sizeof(double), &results, RESULTFIELD*sizeof(double), egress );
      
      /* dump */
      
      for ( j=0;j<datacount;j++ )
	{
	  /*      T     nc-in nc-out nf    mu    V     g     E    */
	  printf("%1.4f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n",
		 p.T,dpts[j*DATAFIELD+NC],
		 results[j*RESULTFIELD+ONC],results[j*RESULTFIELD+ONF], results[j*RESULTFIELD+OMU],
		 results[j*RESULTFIELD+OV],results[j*RESULTFIELD+OG],results[j*RESULTFIELD+OF]  );
	}

      /* reload for next iteration */

      for ( k=0;k<datacount;k++ )
	{
	  if ( results[k*RESULTFIELD+ONF]>0. )
	    {
	      dpts[k*DATAFIELD+MU] = results[k*RESULTFIELD+OMU];
	      dpts[k*DATAFIELD+V] = results[k*RESULTFIELD+OV];
	      dpts[k*DATAFIELD+G] = results[k*RESULTFIELD+OG];
	    }
	}

      /* flush streams */
  
      fflush(stdout);

      /* clear results */

      if (results)
	{
	  free(results);
	  results = NULL;
	}
      
#endif

      /* set temperature for next iteration */

      i++;
      if ( NEAR(LOWT,1.,PREC) )
	{
	  p.T = ((double)i)*INCT;
	}
      else
	{
	  p.T = LOWT + ((double)i)*INCT;
	}

    }
  
  
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
  double vresult;
  double saveg, testg, testj;
  int jflag;
  int rflag = 0;
  char * ecode = NULL;
  int iter = 0;

  p->nc = input[NC];
  p->V = input[V];
  p->mu = input[MU];
  p->g = input[G];


  while ( iter < FAILNUMBER )
    {

      jflag = 0;

      p->mu = mu_finder( mu_wrapper_kondo,p,input[NC] );
      testj = 1./SUMMOR( sce_kondo,p,SUMPREC );
      testg = glam_finder( glam_wrapper_kondo,p,1.0)  ;
      

      if ( NEAR(p->J,testj,PREC) )
	{
	  jflag = 1;
	  if ( NEAR(testg,p->g,PREC) )
	    {
	      break;
	    }
	  else
	    {
	      saveg = p->g;
	    }
	}

      if ( testj > p->J )
	{
	  /* try (0,V) */
	  vresult = bissect_V(p,PREC,p->V);
	}
      else
	{
	  /* try (V,V+inc) */
	  vresult = bissect_V(p,p->V,UPPERLIMITV);
	}
      
      
      /* decision tree */
      
      if ( vresult <= 0. )
	{
	  if ( NEAR(SUMMOR(occupation_kondo_f,p,SUMPREC),1.,PREC*5.)==0 )
	    {
	      p->V = shotgun_shot(0.,2.*input[V]);
	    }
	  if ( jflag )
	    {
	      if ( rflag >= RETRYNUMBER )
		{
		  rflag = 0;
		  p->g = shotgun_shot(LB,UB);		  
		}
	      else
		{
		  rflag++;
		  p->g = saveg;
		}
	    }
	  /*else
	    {
	    p->g = shotgun_shot(LB,UB);
	    }*/
	}
      else
	{
	  p->g = glam_finder(glam_wrapper_kondo,p,1.0)  ;
	  if ( jflag )
	    {
	      if ( rflag >= RETRYNUMBER )
		{
		  rflag = 0;
		  p->V = shotgun_shot(0.,2.*input[V]);
		}
	      else
		{
		  rflag++;
		  if ( fabs(saveg-testg) < fabs(saveg-p->g) )
		    {
		      p->g = saveg;
		      continue;
		    }
		}
	    }
	}

      iter++;

    }
  
  if ( iter>=FAILNUMBER )
    {
      output[ONC] = -1.;
      output[ONF] = -1.;
      output[OMU] = 0.;
      output[OG]  = 0.;
      output[OV] = -1.;
      output[OF] = 0.;
    }
  else
    {
      output[ONC] = SUMMOR( occupation_kondo_c,p,SUMPREC );
      output[ONF] = SUMMOR( occupation_kondo_f,p,SUMPREC );
      output[OMU] = p->mu;
      output[OG]  = p->g;
      output[OV] = p->V;
      /*if I'm using nc/4, then I must mult the energy by 4. But if not, then I should not mult by anything.*/
      output[OF] = 4.*(  SUMMOR(free_energy_kondo,p,SUMPREC) + 3.75*(p->J)*(p->V)*(p->V)  ) - p->g ;
    }

  return (void *)output;
}





int solver_s( struct Parameters * p,
	      double * input,
	      double * output )
{

  char * ecode = NULL;
  int iter = 0;
  double testg, saveg;
  double testj;
  double vresult;
  int jflag;          /* identifies whether initial J value is self-consistent */
  int rflag = 0;      /* identifies the number of times a failed V has resulted in a use of "saveg" */
  
  p->nc = input[NC]; /* fixed. This is the target. */
  p->V = input[V];   /* guess: starting point; should be near the final value */
  p->mu = input[MU]; /* guess: unknown */
  p->g = input[G];   /* guess: starting point */

#ifdef DEBUG
  fprintf(p->fout,"T %f, V %f, mu %f, g %f, nc %f\n",p->T,p->V,p->mu,p->g,p->nc);
  fflush(p->fout);
#endif

  
  while ( iter < FAILNUMBER )
    {

      jflag = 0;
      
#ifdef DEBUG
      fprintf(p->fout,"     ----------------------------------------------------------%d\n     -                         g %1.4e V %1.4e mu %1.4e nc %1.4e\n",iter,p->g,p->V,p->mu,p->nc);
      fflush(p->fout);
#endif


      p->mu = mu_finder( mu_wrapper_kondo,p,input[NC] );
      testj = 1./SUMMOR( sce_kondo,p,SUMPREC );
      testg = glam_finder( glam_wrapper_kondo,p,1.0)  ;
      
#ifdef DEBUG
      fprintf(p->fout,"     -                         mu %1.4e j1 %1.4e\n",p->mu,testj);
      fflush(p->fout);
#endif

      if ( NEAR(p->J,testj,PREC) )
	{
	  jflag = 1;
#ifdef DEBUG
	  fprintf(p->fout,"     -                         J is reached without inc.\n");
	  fflush(p->fout);
#endif
	  if ( NEAR(testg,p->g,PREC) )
	    {
	      break;
	    }
	  else
	    {
	      saveg = p->g;
#ifdef DEBUG
	      fprintf(p->fout,"     -                         g does not match\n");
	      fflush(p->fout);
#endif
	    }
	}

      if ( testj > p->J )
	{
	  /* try (0,V) */
#ifdef DEBUG
	  fprintf(p->fout,"     -                         Bissecting [%1.4e,%1.4e]\n",PREC,p->V);
	  fflush(p->fout);
#endif	  
	  vresult = bissect_V(p,PREC,p->V);
	}
      else
	{
	  /* try (V,V+inc) */
#ifdef DEBUG
	  fprintf(p->fout,"     -                         Bissecting [%1.4e,%1.4e]\n",p->V,UPPERLIMITV);
	  fflush(p->fout);
#endif	  
	  vresult = bissect_V(p,p->V,UPPERLIMITV);
	}
      
      
#ifdef DEBUG
      fprintf(p->fout,"     -                         Vresult %1.4e mu %1.4e\n",vresult,p->mu);
      fflush(p->fout);
#endif
      
      /* decision tree */
      
      if ( vresult <= 0. )
	{
	  if ( NEAR(SUMMOR(occupation_kondo_f,p,SUMPREC),1.,PREC*5.)==0 )
	    {
	      p->V = shotgun_shot(0.,2.*input[V]);
	    }
	  if ( jflag )
	    {
	      if ( rflag >= RETRYNUMBER )
		{
		  rflag = 0;
		  p->g = shotgun_shot(LB,UB);		  
		}
	      else
		{
		  rflag++;
		  p->g = saveg;
		}
	    }
	  /*else
	    {
	      p->g = shotgun_shot(LB,UB);
	      }*/
	}
      else
	{
	  p->g = glam_finder(glam_wrapper_kondo,p,1.0)  ;
#ifdef DEBUG
	  fprintf(p->fout,"     -                         glam finder: %1.4e\n",p->g);
	  fflush(p->fout);
#endif
	  if ( jflag )
	    {
	      if ( rflag >= RETRYNUMBER )
		{
		  rflag = 0;
		  p->V = shotgun_shot(0.,2.*input[V]);
		}
	      else
		{
		  rflag++;
		  if ( fabs(saveg-testg) < fabs(saveg-p->g) )
		    {
		      p->g = saveg;
		      continue;
		    }
		}
	    }
	}
#ifdef DEBUG
      fprintf(p->fout,"     -                         results V %1.4e g %1.4e\n",p->V,p->g);
      fflush(p->fout);
#endif

      iter++;
	  
#ifdef DEBUG
      fprintf(p->fout,"     -                         looping on g\n");
      fflush(p->fout);
#endif
    }

  if ( iter>=FAILNUMBER )
    {
      output[ONC] = -1.;
      output[ONF] = -1.;
      output[OMU] = 0.;
      output[OG]  = 0.;
      output[OV] = -1.;
      output[OF] = 0.;
      return -1;
    }
  else
    {
      output[ONC] = SUMMOR( occupation_kondo_c,p,SUMPREC );
      output[ONF] = SUMMOR( occupation_kondo_f,p,SUMPREC );
      output[OMU] = p->mu;
      output[OG]  = p->g;
      output[OV] = p->V;
      /*if I'm using nc/4, then I must mult the energy by 4. But if not, then I should not mult by anything.*/
      output[OF] = 4.*(  SUMMOR(free_energy_kondo,p,SUMPREC) + 3.75*(p->J)*(p->V)*(p->V)  ) - p->g ;
#ifdef DEBUG
      //fprintf(p->fout,"------- FE: 4 x ( %1.4e + %1.4e ) - %1.4e = %1.4e\n",bande,3.75*(p->J)*(p->V)*(p->V),p->g,output[OF]);
      //fflush(p->fout);
#endif
      return 0;
    }
  
}






double
bissect_V( void * param,
	   double lb,
	   double ub       )
{

  struct Parameters * p = param;
  double lower = lb;
  double upper = ub;
  double res;

  p->V = lower;
  p->mu = mu_finder( mu_wrapper_kondo,p,p->nc );
  res = 1./SUMMOR( sce_kondo,p,SUMPREC );
#ifdef DEBUG
  fprintf(p->fout,"     -                         ...   Lower: %1.4e -> %1.4e\n",p->V,res);
  fflush(p->fout);
#endif	  

  if ( res > p->J )
    return -1.;

  do
    {
      p->V = (lower+upper)/2.;
      p->mu = mu_finder( mu_wrapper_kondo,p,p->nc );
      res = 1./SUMMOR( sce_kondo,p,SUMPREC );
#ifdef DEBUG
      fprintf(p->fout,"     -                         ...   [%1.4e,%1.4e] %1.4e -> %1.4e\n",lower,upper,p->V,res);
      fflush(p->fout);
#endif
      if ( NEAR(res,p->J,PREC) )
	return p->V;      
      if ( res>p->J )
	{
	  upper = p->V;
	}
      else
	{
	  lower = p->V;
	}
    } while ( NEAR(lower,upper,PREC/2.)==0 ) ;

  if ( NEAR(res,p->J,PREC)==0 )
    return -1.*p->V;
  else
    return p->V;
  
}





int glam_looper( double (*gf)(void *, double),
		 double (*mf)(void *, double),
		 void * par,
		 double nimp )
{
  struct Parameters * p = par;
  double runningg[4];
  int g_iter = 0;
  p->g = glam_finder( gf,par,nimp );   /* nf constraint is solved here */
  runningg[1]=-200.;
  runningg[2]=-100.;
  runningg[3]=-50.;
  while ( g_iter < FAILNUMBER )
    {
      runningg[0] = runningg[1];
      runningg[1] = runningg[2];
      runningg[2] = runningg[3];
      p->mu = mu_finder( mf,p,p->nc );
#ifdef DEBUG
      fprintf(p->fout,"     -                                       -> g %1.8e mu %1.8e\n",p->g,p->V);
      fflush(p->fout);
#endif
      p->g = glam_finder( gf,p,nimp );   /* nf constraint is solved here */
      runningg[3] = p->g;	  
      if ( NEAR(runningg[2],p->g,PREC)==1 )
	{
	  break;
	}
      if ( NEAR(runningg[3],runningg[1],PREC)==1 && NEAR(runningg[2],runningg[0],PREC)==1 )
	{
	  p->g = shotgun_shot( LB,UB );
	}
      g_iter++;
    }
  if ( g_iter >= FAILNUMBER )
    {
      return -1;
    }
  return 0;
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
