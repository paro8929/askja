/*********************
ASKJA -- a hybrid monte carlo simulation code 
for SU(N) in arbitrary dimensions

by

Paul Romatschke


v1.0: 2016-12-16

Copyright is granted if you cite the relevant paper(s)

 ********************/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <fstream>
#include <string.h>
#include <structs.h>
#include <complex>
#include <sys/stat.h>
#include <Eigen/Eigenvalues> 
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


//#define filtering 1

using namespace std;
using namespace Eigen;

int VERBOSE =1;

int CONFTAKE=10; //how often to actually take configurations
int NACC=0; //number of accepted configurations
int NAC1=0; //accepted because energy was lowered
int NREC=0; //number of rejected configurations
int TAKEN=0;


int ABORT = 0; //flag signalling to abort asap

int DIM=3; //number of space dimensions, controlled by params file

int NC=2; //number of colors in U(NC)

double BETA=1; //lattice beta

int STEPS=100; //number of steps to take

int UPDATE=10; //store results this often

int UPDATE2=10; //more rare storage occurence

int LEVEL=5; //number of terms to keep in exp(x)=sum_n=0^LEVEL x^n/n!

int RUNCONTINUE=0; //continue run from stored configuration (1) or not (0)

int START=1; //starting time, potentially different if run is being continued

int RANDOM=0; //random initial conditions (1) or fixed seed (0)

double dt=0.1; //time increment in units of lattice spacing

double FILTER=10.; // strength of bandpass filter, standard value 10.

double CRITENERGY; // energy according to which configuration should be kept

complex<double> I(0,1);

data_array xygrid;

SUNmatrix *gens; //u(n) generators, used for initial conditions

sunfield *pi;   //main values
sunfield *u; 

sunfield *Pi;   //updated values
sunfield *U; 


sunfield *uconf; //stored configuration

int Ndust=1;
particle *dust; //particles to measure temperature



void allocate_memory()
{
  printf("===> Info: Allocating memory...\t");
  
  pi = new sunfield[xygrid.SIZE];
  Pi = new sunfield[xygrid.SIZE];
  u = new sunfield[xygrid.SIZE];
  U = new sunfield[xygrid.SIZE];

  uconf = new sunfield[xygrid.SIZE];




  for (long int i=0;i<xygrid.SIZE;i++)
    {
      pi[i].init(xygrid.gridsize.d,NC);
      Pi[i].init(xygrid.gridsize.d,NC);
      u[i].init(xygrid.gridsize.d,NC);
      U[i].init(xygrid.gridsize.d,NC);
      uconf[i].init(xygrid.gridsize.d,NC);
      pi[i].alloc_data();
      Pi[i].alloc_data();
      u[i].alloc_data();
      U[i].alloc_data();
      uconf[i].alloc_data();
    }

  printf(" finished\n");
}


void free_memory()
{
  printf("===> Info: Freeing  memory...\t");
  delete [] pi;
  delete [] Pi;
  delete [] u;
  delete [] U;
  delete [] uconf;

  //delete [] f;
  printf("finished\n");
}


void generators(SUNmatrix* gen)
{
  
    
  
  int num=0;
  gsl_complex dummy1,dummy2,dummy3;
  GSL_SET_COMPLEX(&dummy1, 1, 0); 
  GSL_SET_COMPLEX(&dummy2, 0, 1);
  GSL_SET_COMPLEX(&dummy3, 0, -1);

  //off-diagonal ones
  for (int i=0;i<NC;i++)
    for (int j=i+1;j<NC;j++)
      {
	gen[num].set(i,j,dummy1);
	gen[num].set(j,i,dummy1);
	//printf("setting %i\n",num);
	gen[num+1].set(i,j,dummy3);
	gen[num+1].set(j,i,dummy2);
	//printf("setting %i\n",num+1);
	num+=2;
      }
	 

  
  //diagonal ones
  for (int i=0;i<NC-1;i++)
    {
      gsl_complex dummy1,dummy2;
      GSL_SET_COMPLEX(&dummy1, 1, 0);     
      for (int j=0;j<i+1;j++)
	gen[num].set(j,j,dummy1);

      GSL_SET_COMPLEX(&dummy2, -i-1, 0);
      gen[num].set(i+1,i+1,dummy2);
      num++;
      }

  
  //normalize
  double norm=0;
  for (int num=0;num<NC*NC-1;num++)
    {
      norm=sqrt(0.5/TrR(gen[num]*gen[num]));
      gen[num]*=norm;
    }
}




void initializeU()
{
  printf("===> Info: Initializing gauge fields...\t");
  
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      {
	
	gsl_matrix_complex_set_identity(u[i].M[j].M);
	U[i].M[j]=u[i].M[j];
      }
    
  printf("finished\n");
}

void initializeE(unsigned long int seed)
{
  printf("===> Info: Initializing momenta...\t");
  
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r,seed);

  


  SUNmatrix test(NC);

  //test.print();

  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      {
	
	gsl_matrix_complex_set_zero(test.M);
	
	for (int n=0;n<NC*NC-1;n++)
	  {
	    //printf("i=%li j=%i n=%i rr=%f\n",i,j,n,rannum[n][i*xygrid.gridsize.d+j]);
	    //bandpass(rannum[n][j]);
	    double nummi=gsl_ran_gaussian(r,1.0);
	    //double nummi=gsl_rng_uniform(r);
	    //printf("%i, %g\n",n,nummi);
	    test=test+gens[n]*nummi;
	  }
		
	pi[i].M[j]=test*sqrt(2*NC/BETA);
      }
    
  gsl_rng_free(r);

  printf("finished\n");
}

void loadall(char *outname)
{
  
  char buffer[255];
  sprintf(buffer,"%s.data",outname);

  printf("===> Info: Loading fields...\t");
  
  FILE * pFile;
  pFile = fopen (buffer , "r");
       
  if (pFile == NULL) 
    {
      printf("Error loading fields file %s",buffer);
      printf(" ABORTED!\n");
    }
  else
    {
      for (long int i=0;i<xygrid.SIZE;i++)
	for (int j=0;j<xygrid.gridsize.d;j++)
	  {
	    gsl_matrix_complex_fread(pFile,u[i].M[j].M);
	    gsl_matrix_complex_fread(pFile,U[i].M[j].M);
	    gsl_matrix_complex_fread(pFile,pi[i].M[j].M);
	    gsl_matrix_complex_fread(pFile,Pi[i].M[j].M);

	  }
        printf(" finished successfully\n");
    }
  fclose (pFile);


  
}

void storeconf(char *outname)
{

  printf("==>Info: storing configuration %s\n",outname);
  FILE * pFile;
  pFile = fopen (outname , "w");
  
      for (long int i=0;i<xygrid.SIZE;i++)
	for (int j=0;j<xygrid.gridsize.d;j++)
	  gsl_matrix_complex_fwrite(pFile,uconf[i].M[j].M);
      fclose (pFile);
      

}




void WilsonLoop(double *eigenvalues)
{
  long int pos=0;
  sunfield W;
  W.init(1,NC);
  W.alloc_data();

  W.M[0]=u[pos].M[0];
  for (int i=0;i<xygrid.gridsize.x[0]-1;i++)
    {
      pos=xygrid.getplus(pos,0);
      W.M[0]=W.M[0]*u[pos].M[0];
    }

  /* Eigen test example
     using Eigen::MatrixXd;
  
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = 0.5;

  EigenSolver<MatrixXd> es(m);
  
  complex<double> lambda1,lambda2;
  lambda1= es.eigenvalues()[0];
  lambda2= es.eigenvalues()[1];

  
  printf("l1=%f+i%f\n",real(lambda1),imag(lambda1));
  printf("l2=%f+i%f\n",real(lambda2),imag(lambda2));


   */

  double wilsonexp=0;
  wilsonexp=gsl_complex_abs(Tr(W.M[0]))/NC;
  eigenvalues[0]=wilsonexp;
  /*
  using Eigen::MatrixXcd;
  MatrixXcd m(NC,NC);
  for(int i=0;i<NC;i++)
    for (int j=0;j<NC;j++)
      m(i,j)=GSL_REAL(gsl_matrix_complex_get(W.M[0].M,i,j))+I*GSL_IMAG(gsl_matrix_complex_get(W.M[0].M,i,j));

  ComplexEigenSolver<MatrixXcd> ces;
  ces.compute(m);

  for (int i=0;i<NC;i++)
    {
      eigenvalues[i+1]=arg(ces.eigenvalues()[i]);
      //printf("lambda(%i)=%f+i %f, abs=%f phase=%f\n",i,real(ces.eigenvalues()[i]),imag(ces.eigenvalues()[i]),abs(ces.eigenvalues()[i]),arg(ces.eigenvalues()[i]));
      }*/
  
}

//fourier trafo of magnetic energy
void FourierTrafo(double *coefficients)
{

  long int pos=0;
  int n=xygrid.gridsize.x[0];
  double data[2*n];

  sunfield staple;
  staple.init(1,NC);
  staple.alloc_data();

  double temp=0;

  for (int i=0;i<n;i++)
    {
      gsl_matrix_complex_set_zero (staple.M[0].M);
      temp=0;
      for (int j=0;j<xygrid.gridsize.d;j++)
	for (int k=j+1;k<xygrid.gridsize.d;k++)
	  {
	    long int njplus=xygrid.getplus(i,j);
	    long int nkplus=xygrid.getplus(i,k);
	    staple.M[0]=u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
	    temp+=1-TrR(staple.M[0])/NC;
	  }

      data[2*i]=temp;
      data[2*i+1]=0.0;
      pos=xygrid.getplus(pos,0);
    }
  

  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;
 
  wavetable = gsl_fft_complex_wavetable_alloc (n);
  workspace = gsl_fft_complex_workspace_alloc (n);

  gsl_fft_complex_inverse (data, 1, n, wavetable, workspace);

  for (int i=0;i<n;i++)
    coefficients[i]=sqrt(data[2*i]*data[2*i]+data[2*i+1]*data[2*i+1]);

  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
  
}



double GaussLaw()
{
  double res=0;
  sunfield temp;
  temp.init(1,NC);
  temp.alloc_data();
  
 


  for (long int i=0;i<xygrid.SIZE;i++)
    {
      for (int j=0;j<xygrid.gridsize.d;j++)
	{
	  long int njminus=xygrid.getminus(i,j);
	  temp.M[0]= temp.M[0]+Pi[i].M[j]-transportLeft(Pi[njminus].M[j],u[njminus].M[j]);
	}
     
    }

  res=tr2(temp.M[0]);

  return res;
}

//checking unitarity violation
double UViol()
{
  double res=0;
  sunfield temp,unity;
  temp.init(1,NC);
  temp.alloc_data();
  
  unity.init(1,NC);
  unity.alloc_data();
  
 
  gsl_matrix_complex_set_identity(unity.M[0].M);

  //unity.M[0].print();
  
  

  for (long int i=0;i<xygrid.SIZE;i++)
    {
      for (int j=0;j<xygrid.gridsize.d;j++)
	{	  
	  
	  temp.M[0]=u[i].M[j];
	  temp.M[0]=temp.M[0]*Hermitian(u[i].M[j])-unity.M[0];
	  res+=tr2(temp.M[0]);
	}
     
    }

  return res;
}

double EnergyMagind(int j)
{
  double res=0;
  sunfield temp;
  sunfield staple;

  temp.init(1,NC);
  temp.alloc_data();
  staple.init(1,NC);
  staple.alloc_data();
  

  for (long int i=0;i<xygrid.SIZE;i++)
    {
      //staple.M[0].a=0;
      gsl_matrix_complex_set_zero (staple.M[0].M);
      for (int k=0;k<xygrid.gridsize.d;k++)
	{
	  if (k!=j)
	    {
	      long int njplus=xygrid.getplus(i,j);
	      long int nkplus=xygrid.getplus(i,k);
	      staple.M[0]=u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
	      res+=1-TrR(staple.M[0])/NC;
	    }
	  
	}
    }
  

  return res;
}


double EnergyMag()
{
  double res=0;
  sunfield temp;
  sunfield staple;

  temp.init(1,NC);
  temp.alloc_data();
  staple.init(1,NC);
  staple.alloc_data();
  

  for (long int i=0;i<xygrid.SIZE;i++)
    {
      //staple.M[0].a=0;
      gsl_matrix_complex_set_zero (staple.M[0].M);
      
      for (int j=0;j<xygrid.gridsize.d;j++)
	for (int k=j+1;k<xygrid.gridsize.d;k++)
	  {
	    long int njplus=xygrid.getplus(i,j);
	    long int nkplus=xygrid.getplus(i,k);
	    staple.M[0]=u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
	    res+=1-TrR(staple.M[0])/NC;
	    
	  }
    }
  
  
  return res;
}


double simpleenergyEl(sunfield *myp)
{
  double res=0;
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      res+=tr2(myp[i].M[j]); 

  return res;
}


double EnergyElind(int j)
{
  double res=0;
  
  for (long int i=0;i<xygrid.SIZE;i++)     
    res+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);    
  return res;
}

double EnergyEl()
{
  double res=0;
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      res+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]); 
  return res;
}


//in lattice units
double TotalEnergy()
{
  double res;
  res=2*NC*EnergyMag()+EnergyEl();
  return res;
}


//energy momentum tensor at site i
void Tab(long int i, double *mdata,int Npt)
{
  

  double res=0;
  sunfield temp,Fij,storeme;

  temp.init(1,NC);
  temp.alloc_data();
  Fij.init(xygrid.gridsize.d*xygrid.gridsize.d,NC);
  Fij.alloc_data();

  storeme.init(xygrid.gridsize.d*xygrid.gridsize.d,NC);
  storeme.alloc_data();

  gsl_matrix_complex_set_zero (temp.M[0].M);


  for (int j=0;j<xygrid.gridsize.d;j++)
    for (int k=0;k<xygrid.gridsize.d;k++)
	  {
	    if (k!=j)
	      {
		long int right=xygrid.getplus(i,j);
		long int up=xygrid.getplus(i,k);
		long int left=xygrid.getminus(i,j);
		long int down=xygrid.getminus(i,k);
		long int downright=xygrid.getplus(down,j);
		long int downleft=xygrid.getminus(down,j);
		long int upleft=xygrid.getminus(up,j);
		
		temp.M[0]=u[i].M[j]*u[right].M[k]*Hermitian(u[up].M[j])*Hermitian(u[i].M[k]);
		//if ((j==0)&&(k==1))
		//  mdata[1]+=NC-TrR(temp.M[0]);
		temp.M[0]=temp.M[0]+u[i].M[k]*Hermitian(u[upleft].M[j])*Hermitian(u[left].M[k])*u[left].M[j];
		temp.M[0]=temp.M[0]+Hermitian(u[left].M[j])*Hermitian(u[downleft].M[k])*u[downleft].M[j]*u[down].M[k];
		temp.M[0]=temp.M[0]+Hermitian(u[down].M[k])*u[down].M[j]*u[downright].M[k]*Hermitian(u[i].M[j]);
		storeme.M[j*xygrid.gridsize.d+k]=temp.M[0]/8;	       
	      }
	    
	  }

  for (int j=0;j<xygrid.gridsize.d;j++)
    for (int k=0;k<xygrid.gridsize.d;k++)
      {
	if (k==j)
	  gsl_matrix_complex_set_zero (Fij.M[j*xygrid.gridsize.d+k].M);
	else
	  Fij.M[j*xygrid.gridsize.d+k]=storeme.M[j*xygrid.gridsize.d+k]-storeme.M[k*xygrid.gridsize.d+j];
	//mdata[0]-=tr2(Fij.M[j*xygrid.gridsize.d+k]);
      }
 
  //to be lifted
  for (int jj=2;jj<Npt;jj++)
    mdata[jj]-=0.5*tr2(pi[i].M[jj-2]+Pi[i].M[jj-2]);


  //assume that mdata entries are set to zero initially
  for (int j=0;j<xygrid.gridsize.d;j++)
    {
      //to be lifted
      mdata[0]+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);
      for (int jj=2;jj<Npt;jj++)
      	mdata[jj]+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);

      mdata[1]+=TrI((pi[i].M[j]+Pi[i].M[j])*Fij.M[j*xygrid.gridsize.d+0]);

      //to be lifted
      for (int jj=2;jj<Npt;jj++)
      	mdata[jj]-=2*tr2(Fij.M[j*xygrid.gridsize.d+jj-2]);      




      //printf("testing j=%i f0%i=%g and f1%i=%g\n",j,j,tr2(Fij.M[j*xygrid.gridsize.d+0]),j,tr2(Fij.M[j*xygrid.gridsize.d+1]));
      
      for (int k=j+1;k<xygrid.gridsize.d;k++)
	{
	 

	  //to be lifted

	  mdata[0]-=tr2(Fij.M[j*xygrid.gridsize.d+k]);
	  for (int jj=2;jj<Npt;jj++)
	    mdata[jj]+=tr2(Fij.M[j*xygrid.gridsize.d+k]);
	}
    }

 

}

void updateE(double eps)
{
  sunfield temp;
  sunfield staple;
  int nthreads, tid;
  

#pragma omp parallel shared(xygrid,u,pi,Pi,eps,nthreads) private(temp,staple,tid)
  {
    temp.init(1,NC);
    temp.alloc_data();
    staple.init(1,NC);
    staple.alloc_data();

    tid = omp_get_thread_num();
    if (tid == 0)
      {
	nthreads = omp_get_num_threads();
	//printf("===> Info: Number of threads = %d\n", nthreads);
      }

#pragma omp for schedule(dynamic,1)
    for (long int i=0;i<xygrid.SIZE;i++)
      {      
	for (int j=0;j<xygrid.gridsize.d;j++)
	  {
	    long int njplus=xygrid.getplus(i,j);
	    
	    gsl_matrix_complex_set_zero (staple.M[0].M);
	    
	    for (int k=0;k<xygrid.gridsize.d;k++)
	      {	   
		long int nkplus=xygrid.getplus(i,k);
		long int nkminus=xygrid.getminus(i,k);
		long int njkminus=xygrid.getminus(njplus,k);
		
		
		if (k!=j)
		  {
		    
		    staple.M[0]=staple.M[0]+u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
		    staple.M[0]=staple.M[0]+u[i].M[j]*Hermitian(u[njkminus].M[k])*Hermitian(u[nkminus].M[j])*u[nkminus].M[k];		
		  }
	      }
	    Pi[i].M[j]=pi[i].M[j]-Adjoint(staple.M[0])*eps;
	    
	    
	  }
      }
  }//end of parallel session

}

void takeconfiguration()
{
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      uconf[i].M[j]=u[i].M[j];
  if (VERBOSE)
    printf("===> Info: accepted new configuration, em=%g\n",EnergyMag()*BETA);
}

void loadconfiguration()
{
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      u[i].M[j]=uconf[i].M[j];

  if (VERBOSE)
    printf("===> Info: loaded old configuration, em=%g\n",EnergyMag()*BETA);
}

void loadconf(char *outname)
{
  //printf("===> Info: Loading fields...\t");

  FILE * pFile;
  pFile = fopen (outname , "r");
  
  if (pFile == NULL) 
    {
      printf("Error loading fields file %s",outname);
      printf(" ABORTED!\n");
    }
  else
    {
      for (long int i=0;i<xygrid.SIZE;i++)
	for (int j=0;j<xygrid.gridsize.d;j++)
	  gsl_matrix_complex_fread(pFile,u[i].M[j].M);
      //printf(" finished successfully\n");
    }
  fclose (pFile);

}

void updateU(double eps)
{
  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      U[i].M[j]=mexpi(Pi[i].M[j]*eps,LEVEL)*u[i].M[j];
}

void updateUandcopyDown(double eps)
{
  int nthreads, tid;
#pragma omp parallel shared(xygrid,u,U,Pi,eps,nthreads,LEVEL) private(tid)
  {
#pragma omp for schedule(dynamic,1)
    for (long int i=0;i<xygrid.SIZE;i++)
      {
	for (int j=0;j<xygrid.gridsize.d;j++)
	  {
	    U[i].M[j]=mexpi(Pi[i].M[j]*eps,LEVEL)*u[i].M[j];
	    u[i].M[j]=U[i].M[j];
	    pi[i].M[j]=Pi[i].M[j];
	  }
      }//end of parallel session
  }
}

void copyDown()
{

  for (long int i=0;i<xygrid.SIZE;i++)
    for (int j=0;j<xygrid.gridsize.d;j++)
      {
	u[i].M[j]=U[i].M[j];
	pi[i].M[j]=Pi[i].M[j];
      }
}

void load_parameters(char * paramsfile,struct coord * mygrid)
{

  fstream parameters;
  
 
  parameters.open(paramsfile, ios::in);
  if(parameters.is_open())
    {      
      int dummyint;
      char dummychar[255];
      double dummydouble;
      parameters >> dummychar;
      parameters >> dummychar;

      parameters >> dummyint;
      NC=dummyint;

      parameters >> dummychar;
      parameters >> dummyint;
      DIM=dummyint;
      

      printf("DIM =%i\n",DIM);
      mygrid->d=DIM;
      mygrid->allocate();

      for (int i=0;i<DIM;i++)
	{
	  parameters >> dummychar;
	  parameters >> dummyint;
	  mygrid->x[i]=dummyint;
	}

       parameters >> dummychar;
       parameters >> dummydouble;
       BETA=dummydouble;

       printf("BETA =%f\n",BETA);

       parameters >> dummychar;
       parameters >> dummydouble;
       dt=dummydouble;

       printf("dt =%f\n",dt);

       parameters >> dummychar;
       parameters >> dummyint;
       STEPS=dummyint;

       printf("STEPS =%i\n",STEPS);

       parameters >> dummychar;
       parameters >> dummyint;
       UPDATE=dummyint;

       parameters >> dummychar;
       parameters >> dummyint;
       UPDATE2=dummyint;

       parameters >> dummychar;
       parameters >> dummyint;
       CONFTAKE=dummyint;

       printf("CONFTAKE =%i\n",CONFTAKE);     

       parameters >> dummychar;
       parameters >> dummyint;
       LEVEL=dummyint;

       parameters >> dummychar;
       parameters >> dummyint;
       VERBOSE=dummyint;

       /*
       parameters >> dummychar;
       parameters >> dummydouble;
       ENERGY=dummydouble;*/

       parameters >> dummychar;
       parameters >> dummyint;
       RUNCONTINUE=dummyint;

       parameters >> dummychar;
       parameters >> dummyint;
       RANDOM=dummyint;

       if (RUNCONTINUE)
	 {
	   parameters >> dummychar;
	   parameters >> dummyint;
	   START=dummyint;
	   printf("===> Info: Continuing run at t=%f\n",START*dt);
	   parameters >> dummychar;
	   parameters >> dummyint;
	   TAKEN=dummyint;

	   parameters >> dummychar;
	   parameters >> dummyint;
	   NACC=dummyint;
	   //correct for non-stored parts
	   NACC=(NACC/CONFTAKE)*CONFTAKE;
	   
	   parameters >> dummychar;
	   parameters >> dummyint;
	   NAC1=dummyint;

	   parameters >> dummychar;
	   parameters >> dummyint;
	   NREC=dummyint;

	   printf("INFO: Taken=%i Nacc=%i Nac1=%i Nrec=%i\n",TAKEN,NACC,NAC1,NREC);  
	 }
       else
	 {
	   START=1;
	   TAKEN=0;
	   NACC=0;
	   NAC1=0;
	   NREC=0;	     
	 }

     parameters.close();
    }
  else
    {
      printf("\nERROR: params file %s could not be opened\n",paramsfile);
      ABORT=1;
    }
  
}

void store_parameters(char * infile,char * outfile,int myt, int mytaken)
{

  fstream myin,myout;

  char buffer1[255];
  char buffer2[255];
  sprintf(buffer1,"%s-tmp",infile);
  sprintf(buffer2,"cp %s %s",infile,buffer1);
  system(buffer2);
  
  
  myin.open(buffer1, ios::in);
  myout.open(outfile, ios::out);
  
  if((myin.is_open())&&(myout.is_open()))
    {

      
      int dummyint;
      char dummychar[255];
      double dummydouble;

      myin >> dummychar;
      myout << dummychar << "\n";
      myin >> dummychar;
      myout << dummychar<< "\t\t";

      myin >> dummyint;
      myout << dummyint << "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      
      myin >> dummyint;
      myout << dummyint << "\n";
      
      
      for (int i=0;i<DIM;i++)
	{
	  myin >> dummychar;
	  myout << dummychar<< "\t\t";
	  myin >> dummyint;
	  myout << dummyint<< "\n";
	}
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummydouble;
      myout << dummydouble << "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummydouble;
      myout << dummydouble << "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << STEPS<< "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << UPDATE<< "\n";

      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << UPDATE2<< "\n";

      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << CONFTAKE << "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t";
      myin >> dummyint;
      myout << LEVEL<< "\n";

      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << VERBOSE<< "\n";
           
      myin >> dummychar;
      myout << dummychar<< "\t";
      myin >> dummyint;
      myout << 1<< "\n";
      
      myin >> dummychar;
      myout << dummychar<< "\t\t";
      myin >> dummyint;
      myout << RANDOM<< "\n";
      
      myout << "TIME-STEP\t" << myt << "\n";
      myout << "TAKEN\t\t" << mytaken << "\n";
      
      myout << "NACC\t\t" << NACC << "\n";
      myout << "NAC1\t\t" << NAC1 << "\n";
      myout << "NREC\t\t" << NREC << "\n";

     myin.close();
     myout.close();


     sprintf(buffer2,"rm %s",buffer1);
     system(buffer2);
    }
  else
    {
      printf("\nERROR: one of the params files %s,%s could not be opened\n",infile,outfile);
      ABORT=1;
    }
  
}

int main(int argc, char* argv[])
{


  
  char *paramsfile;
  struct coord mygrid;


  if (argc>1)
    {
      paramsfile=argv[1];
    }
  else
    {
      printf("ERROR please specify a params file name");      
    }

  load_parameters(paramsfile, &mygrid);
  
  

  if (!ABORT)
    {
      printf("SU(%i) simulation in %i space dimensions\n",NC,DIM);

      printf("===> Info: Initializing SU(%i) generators....\t",NC);
      gens=new SUNmatrix[NC*NC-1]; 
      for (int i=0;i<NC*NC-1;i++)
	gens[i].SUNmatrix_init(NC);

      generators(gens);
      

      printf(" finished\n");
      xygrid.init(mygrid);
      xygrid.gridsize.print();
      
      allocate_memory();
      
      
      double em;
      fstream constraints,Wilson;
      
      char buffer[255];
      char outdir[255];
      
      sprintf(outdir,"output-DIM%i",DIM);
      
      
      struct stat st;

      if (stat(outdir,&st) != 0)
	{
	  printf("===> Info: directory %s not present, creating it\n",outdir);
	  mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
      

      char tmpchar[1255];
      double evs[NC+1]; //eigenvalues of Wilson loop plus expectation value
       
      double coefficients[xygrid.gridsize.x[0]]; //Fourier coefficients;


      sprintf(buffer,"%s/constraints-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);

      if (!RUNCONTINUE)
	constraints.open(buffer,ios::out);
      else
	constraints.open(buffer,ios::out | ios::app);
      
      sprintf(buffer,"%s/Wilson-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);
      if (!RUNCONTINUE)
	Wilson.open(buffer,ios::out);
      else
	Wilson.open(buffer,ios::out | ios::app);

      if (!RUNCONTINUE)
	{
	  sprintf(buffer,"#Time\t Gauss-law\t Total Energy\t E_L\t E_M\t Unitary-Violation\n");
	  constraints.write(buffer,strlen(buffer));
	  sprintf(buffer,"#Time\t <W>\t Wilson-Line-Phases\n");
	  Wilson.write(buffer,strlen(buffer));
	}
      
      unsigned long int seed;

 
      if (!RUNCONTINUE)
	{
	  if (RANDOM)
	    seed=time(NULL);
	  else
	    seed=120351616;
	  
	  initializeU();
	  takeconfiguration();
	}
      else
	{

	  seed=time(NULL);
	  sprintf(buffer,"%s/U-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-%i.conf",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,TAKEN-1);
	  printf("===> Info: Attempting to load configuration %s\n",buffer);
	  loadconf(buffer);
	  for (long int i=0;i<xygrid.SIZE;i++)
	    for (int j=0;j<xygrid.gridsize.d;j++)
	      {		
		U[i].M[j]=u[i].M[j];
	      }
	  takeconfiguration();
	}

            
      initializeE(seed);
      CRITENERGY=(EnergyMag()+simpleenergyEl(pi)/(2*NC))*BETA;
      printf("\t INFO: Initial energy =%f\n",CRITENERGY);
      
      updateE(dt/2); //first half step to initialize electric fields.     
      copyDown();
      
      em=EnergyMag();              
      double ee=EnergyEl();
      double gauss=GaussLaw();
      double energy=TotalEnergy();
      WilsonLoop(evs);

      if (!RUNCONTINUE)
	{
	  if (VERBOSE)
	    printf("t=%f Gauss Law =%g Energy=%g Wx=%g EL=%g EM=%g\n",0., gauss,energy/(2*NC)*BETA,evs[0],ee/(2*NC)*BETA,em*BETA);
      
	  sprintf(buffer,"%f\t %g\t %g\t %g\t %f\t %g\n",0.,gauss,energy/(2*NC)*BETA,ee/(2*NC)*BETA,em*BETA,UViol());
	  constraints.write(buffer,strlen(buffer));
	  constraints.flush();
	  sprintf(buffer,"%f\t",0.);
	  Wilson.write(buffer,strlen(buffer));
	  for (int i=0;i<NC+1;i++)
	    {
	      sprintf(buffer,"%g\t",evs[i]);
	      Wilson.write(buffer,strlen(buffer));
	    }
	  sprintf(buffer,"\n");
	  Wilson.write(buffer,strlen(buffer));
	  Wilson.flush();
	}
      
      int timekeep;
      
      
      //dust[1].x[0]=8;
	  
      for (timekeep=START+1;timekeep<STEPS+START+1;timekeep++)
	{
	  

	  updateUandcopyDown(dt);
	  if (!(timekeep%UPDATE2==0))
	    updateE(dt);
	  
	  if (timekeep%UPDATE==0)
	    {
	      double gauss=GaussLaw();
	      double ee=EnergyEl();
	      double em=EnergyMag();
	      double energy=TotalEnergy();
	      double uv=UViol();
	      WilsonLoop(evs);
	      if (VERBOSE)
		printf("t=%f Gauss Law =%g Energy=%g Wx=%g EL=%g EM=%g U-Viol=%g \n",(timekeep)*dt, gauss,energy*BETA/(2*NC),evs[0],ee*BETA/(2*NC),em*BETA,uv);
	      
	      
	      sprintf(buffer,"%f\t %g\t %g\t %g\t %f\t %g\n",timekeep*dt,gauss,energy*BETA/(2*NC),ee*BETA/(2*NC),em*BETA,uv);

	      constraints.write(buffer,strlen(buffer));
	      constraints.flush();
	      
	      sprintf(buffer,"%f\t",timekeep*dt);
	      Wilson.write(buffer,strlen(buffer));	      
	      for (int i=0;i<NC+1;i++)
		{
		  sprintf(buffer,"%g\t",evs[i]);
		  Wilson.write(buffer,strlen(buffer));
		}
	      sprintf(buffer,"\n");
	      Wilson.write(buffer,strlen(buffer));
	      Wilson.flush();
	    }
	  
	  if (timekeep%UPDATE2==0) //throw momenta again
	    {
	      updateE(dt/2); //last half step to synchronize electric fields.     
	      //copyDown();
	      
	      double myenergy=(EnergyMag()+simpleenergyEl(Pi)/(2*NC))*BETA;
	      printf("\t INFO: Final energy =%f\n",myenergy);
	      if (myenergy<CRITENERGY) //accept configuration
		{
		  printf("\t ACCEPT case 1: %f<%f\n",myenergy,CRITENERGY);
		  NACC++;
		  NAC1++;

		  takeconfiguration();

		  if ((NACC%CONFTAKE==0)&&(NACC>0)) //actually store configuration somewhere
		    {	
		      sprintf(buffer,"%s/U-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-%i.conf",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,TAKEN);
		      storeconf(buffer);
		      TAKEN++;
		      
		      sprintf(buffer,"%s/master-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.params",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);     
		      store_parameters(paramsfile,buffer,timekeep,TAKEN);
		      
		      
		    }

		  initializeE(time(NULL));
		  CRITENERGY=(EnergyMag()+simpleenergyEl(pi)/(2*NC))*BETA;
		  printf("\t INFO: initial energy =%f\n",CRITENERGY);
		  
		  updateE(dt/2);
		  //copyDown();
		}
	      else
		{
		  const gsl_rng_type * T;
		  gsl_rng * r;
		  
		  gsl_rng_env_setup();
		  
		  T = gsl_rng_default;
		  r = gsl_rng_alloc (T);
		  
		  gsl_rng_set(r,time(NULL));
		  double nummi=gsl_rng_uniform(r);
		  gsl_rng_free(r);
		  if (nummi<exp(-myenergy+CRITENERGY))//accept
		    {
		      printf("\t ACCEPT case 2: %f<%f\n",nummi,exp(-myenergy+CRITENERGY));
		      NACC++;

		      takeconfiguration();

		      if (NACC%CONFTAKE==0) //actually store configuration somewher
			{
			  sprintf(buffer,"%s/U-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-%i.conf",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,TAKEN);
			  storeconf(buffer);
			  TAKEN++;

			  sprintf(buffer,"%s/master-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.params",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);     
			  store_parameters(paramsfile,buffer,timekeep,TAKEN);
			  
			  
      
			  
			}
		      
		      initializeE(time(NULL));
		      CRITENERGY=(EnergyMag()+simpleenergyEl(pi)/(2*NC))*BETA;
		      printf("\t INFO: initial energy =%f\n",CRITENERGY);

		      
		      updateE(dt/2);
		      //copyDown();
		      
		    }
		  else //reject
		    {
		      printf("\t REJECT %f>%f\n",nummi,exp(-myenergy+CRITENERGY));
		      NREC++;


		      loadconfiguration();
		      initializeE(time(NULL));
		      CRITENERGY=(EnergyMag()+simpleenergyEl(pi)/(2*NC))*BETA;
		      printf("\t INFO: initial energy =%f\n",CRITENERGY);

		      
		      
		      updateE(dt/2);
		      //copyDown();
		      
		    }
		}
	      
	      
	      //sprintf(buffer,"%s/stored-N%i-V%.3f-E%.3f-AT%.4f-DT%.4f-%i",outdir,NC,mygrid.x[0]*AT,ENERGY,AT,dt,timekeep);
	      //storeall(timekeep,buffer,paramsfile);
	      
	    }
	}
      
	
      
      constraints.close();
      Wilson.close();
      
      /*
      sprintf(buffer,"%s/master-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.params",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);     
      store_parameters(paramsfile,buffer,timekeep,TAKEN);
      sprintf(buffer,"%s/master-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i.info",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE);
      constraints.open(buffer,ios::out);
      sprintf(buffer,"#TAKEN\t TRAJ\t CONFTAKE \t acc1/acc2\t acc/tot\t rec/tot\n ");
      constraints.write(buffer,strlen(buffer));
      sprintf(buffer,"%i\t %i\t\t %i\t\t %f\t %f\t %f\n ",TAKEN,NACC+NREC,CONFTAKE,NAC1/(1.0*NACC-NAC1),NACC/(1.0*NACC+NREC),NREC/(1.0*NACC+NREC));
      constraints.write(buffer,strlen(buffer));
      constraints.close();
      */
      //store configurations for later use
      //sprintf(buffer,"%s/stored-N%i-V%.3f-E%.3f-AT%.4f-DT%.4f-%i",outdir,NC,mygrid.x[0]*AT,ENERGY,AT,dt,timekeep);
      
      
      //storeall(timekeep,buffer,paramsfile);
      
    }
  free_memory();
  delete [] gens;

  printf("===>Info: Statistics: acc=%i (%i) rec=%i, rate=%f\n",NACC,NAC1,NREC,(NACC*1.0)/(NACC+NREC));
  printf("===>Took %i configurations\n",TAKEN);
  
  


  if (!ABORT)
    printf("Finished successfully\n");
  else
    printf("Finished with errors\n");
}
