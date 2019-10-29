/**********************************************************************************
 *  sm.c -- Parallel-processing library for multicore cpus. See also sm.h.
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
#include "sm.h"


#include <stdio.h>
#include <signal.h>


/* this is each unit of the shared memory.
 * a child process will have its "next" value set to a position in the shared memory,
 * and when calculation is finished, the child will simply assert record[next]->ready=1; record[next]->value=calc'd value
 */
struct pair {
  int ready;
  double value;
};

struct pairv {
  int * ready;
  char * value;
};

/* the linux-style linked list data struct.
 * its body is a pointer to a unit in the shared memory,
 * as such the linked list provides a way to check the "ready" value.
 * there is no need to store a record unit's index.
 */
struct tasklist {
  struct list_head l;
  struct pair * shdp;
};

struct tasklistv {
  struct list_head l;
  int * ready;
  char * value;
};


int d_interrupt;               /* signal-catching flag, used in dispatchv2 */
int smid;                      /* shared memory ID */
int key=0;                     /* shared memory key */
struct pair * record = NULL;   /* for dispatch */
char * recordv = NULL;         /* for dispatchv's */



/* signal handler, used in dispatchv2 */
void d_ctrl(int i)
{
  d_interrupt = 1;
}




int dispatch( double (*f)(void *,double,double,double), void * params, double * data, int n, double ** result )
{

  struct shmid_ds buf;
  struct tasklist tasks;
  struct tasklist * tmpl;
  struct list_head * lpos, *lbuf;
  int next;
  int count;
  int pid = -1;
  int ret;
  int i;
  
  /* set up results memory */
  *result = malloc(n*sizeof(double));
  if ( *result==NULL )
    return -1;
  bzero(*result,n*sizeof(double));
  
  /* set up shared memory */
  key = getpid();
  smid = shmget(key,n*sizeof(struct pair), 0666 | IPC_CREAT | IPC_EXCL);
  if (smid<0)
    return -2;

  /* attach and clear the shared memory portion */
  record = (struct pair *)shmat(smid,NULL,0);
  if ( record==NULL || record<0 )
    return -3;
  bzero(record,sizeof(struct pair)*n);
  
  /* set up the linked list */
  INIT_LIST_HEAD(&tasks.l);

  /* begin dispatch loop */

  next = -1;
  count = 0;
  for(;;)
    {

      /* check the list */
            
      list_for_each_safe(lpos,lbuf,&tasks.l)
	{
	  tmpl = list_entry(lpos,struct tasklist,l);
	  if ( tmpl->shdp->ready > 0 )
	    {
	      waitpid(tmpl->shdp->ready,NULL,0);
	      list_del(lpos);
	      free(tmpl);
	      count--;
	    }
	}
      
      /* add next task */
      
      if ( next < n-1 )
	{
	  if ( count < MAXTASK )
	    {
	      next++;
	      tmpl = malloc(sizeof(struct tasklist));
	      tmpl->shdp = &(record[next]);
	      list_add(&(tmpl->l),&tasks.l);
	      
	      /* fork */
	      
	      pid = fork();
	      if ( pid < 0 )
		{
		  if ( record )
		    {
		      ret = shmdt( record );
		    }
		  if ( smid>=0 )
		    {
		      ret = shmctl( smid, IPC_RMID, &buf );
		    }
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklist,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  return -10;
		  
		}
	      if ( pid == 0 )
		{
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklist,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  record[next].value = f(params,data[3*next],data[3*next+1],data[3*next+2]);
		  record[next].ready = getpid();
		  return 1;
		}
	      else
		{
		  count++;
		  continue;
		}
	    }
	  else
	    {
	    }
	}
      else
	{
	  if ( count==0 )
	    break;
	}
      usleep(SLEEPTIME);
    }

  /* dump data */

  for (i=0;i<n;i++)
    {
      (*result)[i] = record[i].value;
    }
  
  /* clean up*/

  if ( record )
    {
      ret = shmdt( record );
    }
  if ( smid>=0 )
    {
      ret = shmctl( smid, IPC_RMID, &buf );
    }
  list_for_each_safe(lpos,lbuf,&tasks.l)
    {
      tmpl = list_entry(lpos,struct tasklist,l);
      list_del(lpos);
      free(tmpl);
    }
  
  return 0;

}




/* I'm not sure that this version of the above can be called 'complete'. Testing was done with the above, not this one.
   This version is meant, rather than to simply return a 0 to 'main' so that main can exit cleanly & quickly for a child,
   to use a callback or hook function from main. The callback would perform any necessary cleanup and then exit.
   As such, the callback might also be the hook used for a signal-handler in main.
   However, at the time of writing kpms_rhs, it was far simpler to use the above version of dispatch.
*/

void
dispatchv( void * (*f)(void *,void*), void * params, void * dataptr, int n, int sz, void * resultptr, void (*callback)(int,...))
{

  struct shmid_ds buf;
  struct tasklistv tasks;
  struct tasklistv * tmpl;
  struct list_head * lpos, *lbuf;
  int next=0;
  int count=0;
  int pid = -1;
  int ret;
  int i;
  int recordsz;
  char ** result = (char **)resultptr;
  char * data = (char *)dataptr;

  if (callback==NULL)
    return;
  
    /* set up results memory */
  *result = (char *)malloc(n*sz);
  if ( *result==NULL )
    callback(2,-1,dataptr);
  bzero(*result,n*sz);
  
  /* set up shared memory */
  key = getpid();
  recordsz = sizeof(int)+sz;
  smid = shmget(key,n*recordsz, 0666 | IPC_CREAT | IPC_EXCL);
  if (smid<0)
    callback(3,-2,*result,dataptr);
    
  /* attach and clear the shared memory portion */
  recordv = (char *)shmat(smid,NULL,0);
  if ( recordv==NULL || recordv<0 )
    callback(3,-3,*result,dataptr);
  bzero(recordv,recordsz*n);

  /* set up the linked list */
  INIT_LIST_HEAD(&tasks.l);

  /* begin dispatch loop */
  
  next = -1;
  count = 0;
  for(;;)
    {

      /* check the list */
            
      list_for_each_safe(lpos,lbuf,&tasks.l)
	{
	  tmpl = list_entry(lpos,struct tasklistv,l);
	  if ( (*(tmpl->ready)) > 0 )
	    {
	      waitpid(*(tmpl->ready),NULL,0);
	      list_del(lpos);
	      free(tmpl);
	      count--;
	    }
	}
      
      /* add next task */
      
      if ( next < n-1 )
	{
	  if ( count < MAXTASK )
	    {
	      next++;
	      tmpl = malloc(sizeof(struct tasklistv));
	      /*tmpl->shdp = &(recordv[next]);*/
	      tmpl->ready = (int *)(recordv+next*recordsz);
	      tmpl->value = (char *)(recordv+next*recordsz+sizeof(int));
	      list_add(&(tmpl->l),&tasks.l);
	      
	      /* fork */
	      
	      pid = fork();
	      if ( pid < 0 )
		{
		  if ( recordv )
		    {
		      ret = shmdt( recordv );
		    }
		  if ( smid>=0 )
		    {
		      ret = shmctl( smid, IPC_RMID, &buf );
		    }
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklistv,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  callback(3,-10,*result,dataptr);
		  
		}
	      if ( pid == 0 )
		{
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklistv,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  /* i assume here that the results are written back into data contiguously and that the pointer returned by f is to this area of data*/
		  bcopy((char *)f(params,&data[next*sz]),(recordv+next*recordsz+sizeof(int)),sz);
		  *((int *)(recordv+next*recordsz)) = getpid();
		  callback(3,1,dataptr,*result);
		}
	      else
		{
		  count++;
		  continue;
		}
	    }
	  else
	    {
	      /* is there something I wanted to do here ??? */
	    }
	}
      else
	{
	  if ( count==0 )
	    break;
	}
      usleep(SLEEPTIME);
    }

  /* dump data */

  for (i=0;i<n;i++)
    {
      bcopy( (char *)(recordv+i*recordsz+sizeof(int)) , (char *)((*result)+i*sz) , sz );
    }

  /* clean up*/

  if ( record )
    {
      ret = shmdt( recordv );
    }
  if ( smid>=0 )
    {
      ret = shmctl( smid, IPC_RMID, &buf );
    }
  list_for_each_safe(lpos,lbuf,&tasks.l)
    {
      tmpl = list_entry(lpos,struct tasklistv,l);
      list_del(lpos);
      free(tmpl);
    }

  callback(1,0);

}




/* MUST UNREGISTER THE SIGNAL HANDLER BEFORE CALLING THE EXIT-HOOK */
  
void
dispatchv2( struct d_Function dfunc, int n, void * dataptr, int dsz, void * resultptr, int rsz, void (*callback)(int,...))
{

  struct shmid_ds buf;
  struct tasklistv tasks;
  struct tasklistv * tmpl;
  struct list_head * lpos, *lbuf;
  int next=0;
  int count=0;
  int pid = -1;
  int ret;
  int i;
  int recordsz;
  char ** result = (char **)resultptr;
  char * data = (char *)dataptr;
  char * func_return;

  if (callback==NULL)
    return;

  struct sigaction action, oaction;
  sigfillset( &action.sa_mask );
  action.sa_handler = d_ctrl;
  action.sa_flags = 0;
  d_interrupt = 0;
  if ( sigaction(SIGQUIT,&action,&oaction)==-1 )
    {
      /* something here for fail */
    }
  
    /* set up results memory */
  *result = (char *)malloc(n*rsz);
  if ( *result==NULL )
    callback(2,-1,dataptr);
  bzero(*result,n*rsz);
  
  /* set up shared memory */
  key = getpid();
  recordsz = sizeof(int)+rsz;
  smid = shmget(key,n*recordsz, 0666 | IPC_CREAT | IPC_EXCL);
  if (smid<0)
    callback(3,-2,*result,dataptr);
    
  /* attach and clear the shared memory portion */
  recordv = (char *)shmat(smid,NULL,0);
  if ( recordv==NULL || recordv<0 )
    callback(3,-3,*result,dataptr);
  bzero(recordv,recordsz*n);

  func_return = malloc(rsz);
  if ( func_return==NULL )
    callback(3,-4,*result,dataptr);
  bzero(func_return,rsz);

  /* set up the linked list */
  INIT_LIST_HEAD(&tasks.l);

  /* begin dispatch loop */
  
  next = -1;
  count = 0;
  for(;;)
    {

      /* check the list */
            
      list_for_each_safe(lpos,lbuf,&tasks.l)
	{
	  tmpl = list_entry(lpos,struct tasklistv,l);
	  if ( (*(tmpl->ready)) > 0 )
	    {
	      waitpid(*(tmpl->ready),NULL,0);
	      list_del(lpos);
	      free(tmpl);
	      count--;
	    }
	}

      if ( d_interrupt )
	{
	  if ( record )
	    {
	      ret = shmdt( recordv );
	    }
	  if ( smid>=0 )
	    {
	      ret = shmctl( smid, IPC_RMID, &buf );
	    }
	  list_for_each_safe(lpos,lbuf,&tasks.l)
	    {
	      tmpl = list_entry(lpos,struct tasklistv,l);
	      list_del(lpos);
	      free(tmpl);
	    }

	  callback(3,-11,*result,dataptr);

	}
      /* add next task */

      if ( next < n-1 )
	{
	  if ( count < MAXTASK )
	    {
	      next++;
	      tmpl = malloc(sizeof(struct tasklistv));
	      /*tmpl->shdp = &(recordv[next]);*/
	      tmpl->ready = (int *)(recordv+next*recordsz);
	      tmpl->value = (char *)(recordv+next*recordsz+sizeof(int));
	      list_add(&(tmpl->l),&tasks.l);
	      
	      /* fork */
	      
	      pid = fork();
	      if ( pid < 0 )
		{
		  if ( recordv )
		    {
		      ret = shmdt( recordv );
		    }
		  if ( smid>=0 )
		    {
		      ret = shmctl( smid, IPC_RMID, &buf );
		    }
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklistv,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  callback(3,-10,*result,dataptr);
		  
		}
	      if ( pid == 0 )
		{
		  list_for_each_safe(lpos,lbuf,&tasks.l)
		    {
		      tmpl = list_entry(lpos,struct tasklistv,l);
		      list_del(lpos);
		      free(tmpl);
		    }
		  /* i assume here that the results are written back into data contiguously and that the pointer returned by f is to this area of data*/
		  bcopy((char *)(dfunc.f)(dfunc.param,&data[next*dsz],func_return),(recordv+next*recordsz+sizeof(int)),rsz);
		  *((int *)(recordv+next*recordsz)) = getpid();
		  free(func_return);
		  callback(3,1,dataptr,*result);
		}
	      else
		{
		  count++;
		  continue;
		}
	    }
	  else
	    {
	      /* what did you want to do here?? */
	    }
	}
      else
	{
	  if ( count==0 )
	    break;
	}
      usleep(SLEEPTIME);
    }

  /* dump data */

  for (i=0;i<n;i++)
    {
      bcopy( (char *)(recordv+i*recordsz+sizeof(int)) , (char *)((*result)+i*rsz) , rsz );
    }

  /* clean up*/

  if ( record )
    {
      ret = shmdt( recordv );
    }
  if ( smid>=0 )
    {
      ret = shmctl( smid, IPC_RMID, &buf );
    }
  list_for_each_safe(lpos,lbuf,&tasks.l)
    {
      tmpl = list_entry(lpos,struct tasklistv,l);
      list_del(lpos);
      free(tmpl);
    }

  callback(1,0);

}
