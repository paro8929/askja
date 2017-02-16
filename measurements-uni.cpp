#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <fstream>
#include <string.h>
#include <structs.h>
#include <complex>
#include <vector> 
#include <sys/stat.h>
#include <Eigen/Eigenvalues> 
#include <gsl/gsl_fft_complex.h>

using namespace std;
using namespace Eigen;

int ABORT = 0; //flag signalling to abort asap

int VERBOSE =0;

int GET_EV_DIST=0;
int GET_P5_CORR=0;
int GET_ORDERED_WILSON=0;

int CONFTAKE=10; //how often to actually take configurations
int NACC=0; //number of accepted configurations
int NAC1=0; //accepted because energy was lowered
int NREC=0; //number of rejected configurations
int TAKEN=0;

int SKIP=0;
int STRIDE=1;

int largeD=0;

int DIM=3; //number of space dimensions, controlled by params file

int NC=2; //number of colors in SU(NC)

double BETA=1.0; //lattice beta

int STEPS=100; //number of steps to take

int UPDATE=10; //store results this often

int UPDATE2=10; //more rare storage occurence

int LEVEL=5; //number of terms to keep in exp(x)=sum_n=0^LEVEL x^n/n!

int RUNCONTINUE=0; //continue run from stored configuration (1) or not (0)

int START=1; //starting time, potentially different if run is being continued

int RANDOM=0; //random initial conditions (1) or fixed seed (0)

double dt=0.1; //time increment in units of lattice spacing

complex<double> I(0,1);

data_array xygrid;

SUNmatrix *gens; //su(n) generators, used for initial conditions

sunfield *pi;   //main values
sunfield *u; 

sunfield *Pi;   //updated values
sunfield *U; 



void allocate_memory()
{
  printf("===> Info: Allocating memory...\t");
  
  pi = new sunfield[xygrid.SIZE];
  Pi = new sunfield[xygrid.SIZE];
  u = new sunfield[xygrid.SIZE];
  U = new sunfield[xygrid.SIZE];

   for (long int i=0;i<xygrid.SIZE;i++)
    {
      pi[i].init(xygrid.gridsize.d,NC);
      Pi[i].init(xygrid.gridsize.d,NC);
      u[i].init(xygrid.gridsize.d,NC);
      U[i].init(xygrid.gridsize.d,NC);
      pi[i].alloc_data();
      Pi[i].alloc_data();
      u[i].alloc_data();
      U[i].alloc_data();
    }

  printf(" finished\n");
}

static bool sort_using_greater_than(double u, double v)
{
   return u > v;
}


void free_memory()
{
  printf("===> Info: Freeing  memory...\t");
  delete [] pi;
  delete [] Pi;
  delete [] u;
  delete [] U;

  //delete [] f;
  printf("finished\n");
}


void generators(SUNmatrix* gen)
{
  //gen = new SUNmatrix[NC*NC-1];  //to be done outside
  //first diagonal generators;
  //for (int i=0;i<NC*NC-1;i++)
  // gen[i].SUNmatrix_init(NC);//to be done outside
    
  
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


void getevs(double *w,double *p)
{
  long int pos=0;
  long int xpos=0;
  sunfield W;
  W.init(1,NC);
  W.alloc_data();
  
  struct coord co;
  co.d=DIM;
  co.allocate();
  for (int dd=0;dd<DIM;dd++)
    co.x[dd]=0;

  for (int i=0;i<NC;i++)
    {
      w[i]=0;
      p[i]=0;
    }

  //Wilson loop
  for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0];co.x[0]++)
    {
      std::vector<double> myvector (w, w+NC);
      pos=xygrid.getind(&co);
      W.M[0]=u[pos].M[1];

      for (co.x[1]=0;co.x[1]<xygrid.gridsize.x[1]-1;co.x[1]++)
	{
	  pos=xygrid.getplus(pos,1);
	  W.M[0]=W.M[0]*u[pos].M[1];
	}

      //
      //if (VERBOSE)
      //	printf("x=%i Wx=%f\n",co.x[0],gsl_complex_abs(Tr(W.M[0]))/NC);
      

      using Eigen::MatrixXcd;
      MatrixXcd m(NC,NC);
      for(int i=0;i<NC;i++)
	for (int j=0;j<NC;j++)
	  m(i,j)=GSL_REAL(gsl_matrix_complex_get(W.M[0].M,i,j))+I*GSL_IMAG(gsl_matrix_complex_get(W.M[0].M,i,j));
      
      ComplexEigenSolver<MatrixXcd> ces;
      ces.compute(m);
      double winding=0;
      for (int i=0;i<NC;i++)
	{
	  myvector[i]=arg(ces.eigenvalues()[i]);
	  //printf("x0=%i lambda(%i)=%f+i %f, abs=%f phase=%f\n",co.x[0],i,real(ces.eigenvalues()[i]),imag(ces.eigenvalues()[i]),abs(ces.eigenvalues()[i]),arg(ces.eigenvalues()[i]));
	   winding+=myvector[i]/(2*M_PI);
	}

      if ((winding>0.1)||(winding<-0.1)) //if we have non-trivial winding number, recenter to zero
	for (int i=0;i<NC;i++)
	  myvector[i]-=winding*2*M_PI/NC;

      // using default comparison (operator <):
      std::sort (myvector.begin(), myvector.begin()+NC);
      for (int i=0;i<NC;i++)
	{
	  w[i]+=myvector[i];
	  //	printf("i=%i w[i]=%f v=%f\n",i,w[i],myvector[i]);
	}
      
     
    }
  for (int dd=0;dd<DIM;dd++)
    co.x[dd]=0;

  
  //Polyakov loop
  for (co.x[1]=0;co.x[1]<xygrid.gridsize.x[1];co.x[1]++)
    {
      
      pos=xygrid.getind(&co);
      W.M[0]=u[pos].M[0];
      std::vector<double> myvector (p, p+NC);
      for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0]-1;co.x[0]++)
	{
	  pos=xygrid.getplus(pos,0);
	  W.M[0]=W.M[0]*u[pos].M[0];
	}

      using Eigen::MatrixXcd;
      MatrixXcd m(NC,NC);
      for(int i=0;i<NC;i++)
	for (int j=0;j<NC;j++)
	  m(i,j)=GSL_REAL(gsl_matrix_complex_get(W.M[0].M,i,j))+I*GSL_IMAG(gsl_matrix_complex_get(W.M[0].M,i,j));
      
      ComplexEigenSolver<MatrixXcd> ces;
      ces.compute(m);
      double winding=0;
  

      
      for (int i=0;i<NC;i++)
	{
	  myvector[i]=arg(ces.eigenvalues()[i]);
	  //printf("x0=%i lambda(%i)=%f+i %f, abs=%f phase=%f\n",co.x[0],i,real(ces.eigenvalues()[i]),imag(ces.eigenvalues()[i]),abs(ces.eigenvalues()[i]),arg(ces.eigenvalues()[i]));
	  winding+=myvector[i]/(2*M_PI);
	}

      if ((winding>0.1)||(winding<-0.1)) //if we have non-trivial winding number, recenter to zero
	for (int i=0;i<NC;i++)
	  myvector[i]-=winding*2*M_PI/NC;

      // using default comparison (operator <):
      std::sort (myvector.begin(), myvector.begin()+NC);
      for (int i=0;i<NC;i++)
	{
	  p[i]+=myvector[i];
	  //	printf("i=%i w[i]=%f v=%f\n",i,w[i],myvector[i]);
	}
    }
  

  for (int i=0;i<NC;i++)
    {
      w[i]/=xygrid.gridsize.x[0];
      p[i]/=xygrid.gridsize.x[1];
    }
  
}

void WilsonLoop(double *eigenvalues,int dir)
{


  long int pos=0;
  long int xpos=0;
  sunfield W;
  W.init(1,NC);
  W.alloc_data();

  struct coord co;
  co.d=DIM;
  co.allocate();
  for (int dd=0;dd<DIM;dd++)
    co.x[dd]=0;

  double wilsonexp=0;


  for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0];co.x[0]++)
    {
      
      pos=xygrid.getind(&co);
      W.M[0]=u[pos].M[dir];

      for (co.x[dir]=0;co.x[dir]<xygrid.gridsize.x[dir]-1;co.x[dir]++)
	{
	  pos=xygrid.getplus(pos,dir);
	  W.M[0]=W.M[0]*u[pos].M[dir];
	}

      //
      //if (VERBOSE)
      //	printf("x=%i Wx=%f\n",co.x[0],gsl_complex_abs(Tr(W.M[0]))/NC);
      
      wilsonexp+=gsl_complex_abs(Tr(W.M[0]))/NC;
     
     
    }

 
  eigenvalues[0]=(wilsonexp/xygrid.gridsize.x[0]);
  
  
}

void scalarLoopcorr(double *eigenvalues)
{


  long int pos=0;
  long int xpos=0;
  sunfield W;

  double mytr[xygrid.gridsize.x[0]];
  //double mytrdag[xygrid.gridsize.x[0]];
  long int slice=xygrid.helparr[0];
  
  W.init(xygrid.gridsize.x[0],NC);
  W.alloc_data();

  
  
  struct coord co;
  co.d=DIM;
  co.allocate();
  for (int dd=0;dd<DIM;dd++)
    co.x[dd]=0;

  double wilsonexp=0;

  
  for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0];co.x[0]++)
    {
      pos=xygrid.getind(&co);
      W.M[co.x[0]]=u[pos].M[DIM-1];
      mytr[co.x[0]]=TrR(W.M[co.x[0]])/NC;
      //mytrdag[co.x[0]]=TrR(Hermitian(W.M[co.x[0]]))/NC;
      //printf("testing %f + i %f \n",Tr(W.M[co.x[0]]).dat[0],TrR(W.M[co.x[0]]));
      //
      for (long int sl=pos+1;sl<pos+slice;sl++)
	mytr[co.x[0]]+=TrR(u[sl].M[DIM-1])/NC;
	  //mytrdag[co.x[0]]+=TrR(Hermitian(u[sl].M[DIM-1]))/NC;
	
      //W.M[co.x[0]]=W.M[co.x[0]]+u[sl].M[DIM-1];

      //printf("P5[%i]=%f\n",co.x[0],gsl_complex_abs(Tr(W.M[co.x[0]]))/(NC*slice));
    }

  for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0];co.x[0]++)
    {
      //eigenvalues[co.x[0]]=gsl_complex_abs(gsl_complex_mul(Tr(Hermitian(W.M[co.x[0]])),Tr(W.M[0])))/(NC*NC*slice*slice);
      eigenvalues[co.x[0]]=mytr[co.x[0]]/slice;
      if (VERBOSE)
	printf("%i\t %f\n",co.x[0],eigenvalues[co.x[0]]);
    }
    
}

void PolyakovLoop(double *eigenvalues)
{
  long int pos=0;
  long int xpos=0;
  sunfield W;
  W.init(1,NC);
  W.alloc_data();

  struct coord co;
  co.d=DIM;
  co.allocate();
  for (int dd=0;dd<DIM;dd++)
    co.x[dd]=0;

  double wilsonexp=0;

  for (co.x[1]=0;co.x[1]<xygrid.gridsize.x[1];co.x[1]++)
    {
      
      pos=xygrid.getind(&co);
      W.M[0]=u[pos].M[0];
      for (co.x[0]=0;co.x[0]<xygrid.gridsize.x[0]-1;co.x[0]++)
	{
	  pos=xygrid.getplus(pos,0);
	  W.M[0]=W.M[0]*u[pos].M[0];
	}

      // if (VERBOSE)
      //	printf("x=%i W0=%f\n",co.x[1],gsl_complex_abs(Tr(W.M[0]))/NC);
      
      wilsonexp+=gsl_complex_abs(Tr(W.M[0]))/NC;
      

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
	}
      */
      
    }

 
  eigenvalues[0]=(wilsonexp/xygrid.gridsize.x[1]);
  
  
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
void Tab(long int i, double *mdata)
{
  mdata[0]=0;
  mdata[1]=0;
  mdata[2]=0;
  mdata[3]=0;

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
 

  mdata[2]-=0.5*tr2(pi[i].M[0]+Pi[i].M[0]);
  mdata[3]-=0.5*tr2(pi[i].M[1]+Pi[i].M[1]);

  //assume that mdata entries are set to zero initially
  for (int j=0;j<xygrid.gridsize.d;j++)
    {
      mdata[0]+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);
      mdata[1]+=TrI((pi[i].M[j]+Pi[i].M[j])*Fij.M[j*xygrid.gridsize.d+0]);
      mdata[2]-=2*tr2(Fij.M[j*xygrid.gridsize.d+0]);
      mdata[3]-=2*tr2(Fij.M[j*xygrid.gridsize.d+1]);
      //printf("testing j=%i f0%i=%g and f1%i=%g\n",j,j,tr2(Fij.M[j*xygrid.gridsize.d+0]),j,tr2(Fij.M[j*xygrid.gridsize.d+1]));
      mdata[2]+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);
      mdata[3]+=0.25*tr2(pi[i].M[j]+Pi[i].M[j]);
      for (int k=j+1;k<xygrid.gridsize.d;k++)
	{
	  mdata[0]-=tr2(Fij.M[j*xygrid.gridsize.d+k]);
	  mdata[2]+=tr2(Fij.M[j*xygrid.gridsize.d+k]);
	  mdata[3]+=tr2(Fij.M[j*xygrid.gridsize.d+k]);
	}
    }

}

void load_measure_params(char * paramsfile)
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

      parameters >> dummychar;
      parameters >> dummyint;
      SKIP=dummyint;

      parameters >> dummychar;
      parameters >> dummyint;
      STRIDE=dummyint;

      parameters >> dummychar;
      parameters >> dummyint;
      GET_EV_DIST=dummyint;

      parameters >> dummychar;
      parameters >> dummyint;
      GET_P5_CORR=dummyint;

      parameters >> dummychar;
      parameters >> dummyint;
      GET_ORDERED_WILSON=dummyint;


      parameters >> dummychar;
      parameters >> dummyint;
      VERBOSE=dummyint;
      
      
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


void PlaqExp(double *temporal,double *spatial)
{
  *(temporal)=0;
  *(spatial)=0;
  
  //double res=0;
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
      
      for (int j=1;j<xygrid.gridsize.d;j++)
	for (int k=j+1;k<xygrid.gridsize.d;k++)
	  {
	    long int njplus=xygrid.getplus(i,j);
	    long int nkplus=xygrid.getplus(i,k);
	    staple.M[0]=u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
	    *(spatial)+=1-TrR(staple.M[0])/NC;
	    
	  }

      int j=0;
      for (int k=j+1;k<xygrid.gridsize.d;k++)
	{
	  long int njplus=xygrid.getplus(i,j);
	  long int nkplus=xygrid.getplus(i,k);
	  staple.M[0]=u[i].M[j]*u[njplus].M[k]*Hermitian(u[nkplus].M[j])*Hermitian(u[i].M[k]);
	  *(temporal)+=1-TrR(staple.M[0])/NC;	  
	}
    }
  
  
  
}


int main(int argc, char* argv[])
{
  
  char *paramsfile;
  char *mparamsfile;
  struct coord mygrid;

  if (argc>2)
    {
      paramsfile=argv[1];
      mparamsfile=argv[2];            
    }
  else
    {
      printf("ERROR please specify BOTH a params file and a measurement params file name");      
    }

  load_parameters(paramsfile, &mygrid);
  load_measure_params(mparamsfile);
  
  
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
      //xygrid.gridsize.print();
      
      allocate_memory();
      
      char buffer[255];
      char outdir[255];
      
      sprintf(outdir,"rav-DIM%i",DIM);

      struct stat st;
      if (stat(outdir,&st) != 0)
	{
	  printf("===> Info: directory %s not present, creating it\n",outdir);
	  mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

      int vol=mygrid.x[0];

      double evs[NC+1],evsW[NC+1]; //eigenvalues of W & P loop plus expectation value
      double P5[xygrid.gridsize.x[0]];
      double Wssum=0;
      
      
      double Psum=0;
      //double W2sum=0;
      //double Zsum=0;
      double Osum[4];
      for (int i=0;i<4;i++)
	Osum[i]=0.0;


	
      double Utsum=0;
      double Ussum=0;
      double P5sum[xygrid.gridsize.x[0]];

      double Wdist[NC];
      double Pdist[NC];

      double width=M_PI/(2*NC); 
      int bins=floor(2*M_PI/width);
      int Whist[bins],Phist[bins];

      
      if (GET_EV_DIST)
	{
	  
	  for (int i=0;i<NC;i++)
	    {
	      Wdist[i]=0;
	      Pdist[i]=0;
	    }

	  for (int i=0;i<bins;i++)
	    {
	      Whist[i]=0;
	      Phist[i]=0;
	    }
	}
      
      //double P5totsum=0;


      for (int x0=0;x0<xygrid.gridsize.x[0];x0++)
	P5sum[x0]=0;
      
      //int TAKEN=12;
      fstream rav;
      sprintf(buffer,"%s/Nrav-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-SKIP%i-STRIDE%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,SKIP,STRIDE);
      rav.open(buffer,ios::out);

      rav << "#conf" << "\t";
      rav << "P0" << "\t";
      rav << "d(P0)" << "\t";
      rav << "P1" << "\t";
      rav << "dP1" << "\t";
      rav << "U0" << "\t";
      rav << "dU0" << "\t";
      rav << "U1" << "\t";
      rav << "dU1" << "\t";
      rav << "N0" << "\t";
      rav << "N1" << "\t";      
      rav << "BETA" << "\t";
      rav << endl;
      
      fstream P5corr,Hist;
      fstream W2;

      if (GET_ORDERED_WILSON)
	{
	  sprintf(buffer,"%s/WO-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-SKIP%i-STRIDE%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,SKIP,STRIDE);
	  W2.open(buffer,ios::out);
	}

      if (GET_P5_CORR)
	{
	  sprintf(buffer,"%s/P5corr-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-SKIP%i-STRIDE%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,SKIP,STRIDE);
	  P5corr.open(buffer,ios::out);
	}

      if (GET_EV_DIST)
	{
	  sprintf(buffer,"%s/Hist-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-SKIP%i-STRIDE%i.dat",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,SKIP,STRIDE);
	  Hist.open(buffer,ios::out);
	}
      

      sprintf(outdir,"output-DIM%i",DIM);

      printf("===> Info: opening files finished\n");


      //int largeD=0;
      for (int dir=0;dir<xygrid.gridsize.d;dir++)
	{
	  if (xygrid.gridsize.x[dir]>1)
	    largeD++;
	}
      //printf("found %i large D\n",largeD);
	      
      int NAV=0;
      
      for (int conf=SKIP;conf<TAKEN;conf+=STRIDE)//first pass
	{
	  sprintf(buffer,"%s/U-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-%i.conf",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,conf);
	  loadconf(buffer);
	  if (largeD>0)
	    PolyakovLoop(evs);
	  if (largeD>1)
	    WilsonLoop(evsW,1);
	  if (GET_P5_CORR)
	    scalarLoopcorr(P5);
	  double Ut,Us;
	  PlaqExp(&Ut,&Us);

	  if (GET_EV_DIST)
	    {
	      getevs(Wdist,Pdist);

	      for (int i=0;i<NC;i++)
		{
		  int gotme=floor((Wdist[i]+M_PI)/width);
		  if (gotme>bins)
		    printf("ERRROR in histogram binning!\n");
		  else
		    Whist[gotme]++;

		  gotme=floor((Pdist[i]+M_PI)/width);
		  if (gotme>bins)
		    printf("ERRROR in histogram binning!\n");
		  else
		    Phist[gotme]++;
		}
	    }

	  //double em=EnergyMag();
	  /*double sumr=0,sumi=0;
	  for (int cnt=0;cnt<NC;cnt++)
	    {
	      sumr+=cos(evs[cnt+1]);
	      sumi+=sin(evs[cnt+1]);
	      }*/
	  if (VERBOSE)
	    printf("conf=%i Pt=%f Wx=%f Ut=%f Us=%f\n",conf,evs[0],evsW[0],Ut*BETA/vol,Us*BETA/vol);
	  //printf("conf=%i Wx=%f Plaq=%f\n",conf,evs[0],2*NC*em/AT);
	  Wssum+=evsW[0];
	  Psum+=evs[0];
	  double Otemp[largeD];
	  for (int i=0;i<largeD;i++)
	    Otemp[i]=0;

	  if (largeD>0)
	    Otemp[0]=evs[0];
	  if (largeD>1)
	    Otemp[1]=evsW[0];

	  //W2sum+=evs[0]*evs[0];
	  //W2sum+=evsW[0]*evsW[0];
	  //if (evs[0]>0.5)
	  //  Zsum+=1.0;
	  //if (evsW[0]>0.5)
	  //  Zsum+=1.0;	  

	  if (largeD>2)
	    {
	      for (int dir=2;dir<largeD;dir++)
		{
		  WilsonLoop(evsW,dir);
		  Otemp[dir]=evsW[0];
		  //W2sum+=evsW[0]*evsW[0];
		  //if (evsW[0]>0.5)
		  //  Zsum+=1.0;
		}
	    }
	
	      
	  std::vector<double> myvector (Otemp, Otemp+largeD);
	  std::sort (myvector.begin(), myvector.begin()+largeD);

	  for (int i=0;i<largeD;i++)
	    {
	      if (VERBOSE)
		printf("i=%i v=%f\n",i,myvector[i]);
	      Osum[i]+=myvector[i];
	    }
	  
	  //Plaqsum+=em*BETA/vol;
	  Utsum+=Ut*BETA/vol;
	  Ussum+=Us*BETA/vol;

	  for (int x0=0;x0<xygrid.gridsize.x[0];x0++)
	    {
	      P5sum[x0]+=P5[x0];
	      //P5totsum+=P5[x0];
	    }
	  

	  NAV++;
	  printf("Conf %i first pass done!\n",conf);
	  
	  if (conf>SKIP)
	    rav << conf << "\t" << Psum/NAV << "\t 0.\t" << Wssum/NAV << "\t 0.\t" << (Utsum)/NAV << "\t 0.\t" << (Ussum)/NAV << "\t 0.\n";
	  //rav << conf << "\t" << Psum/(conf-SKIP+1) << "\t 0.\t" << Wssum/(conf-SKIP+1) << "\t 0.\t" << Plaqsum/(conf-SKIP+1) << "\t 0.\n";
	  //printf("plaqsum=%f\n",Plaqsum);
	}

      /*
      if (GET_EV_DIST)
	{
	  for (int i=0;i<bins;i++)
	    {
	      //Whist[i]/=(TAKEN-SKIP);
	      printf("i=%i w=%f\n",i,Whist[i]/(1.0*TAKEN-SKIP));
	    }
	    }*/
      
      
      Psum/=NAV;
      Wssum/=NAV;
      Utsum/=NAV;
      Ussum/=NAV;
      for (int x0=0;x0<xygrid.gridsize.x[0];x0++)
	P5sum[x0]/=NAV;

      for (int i=0;i<largeD;i++)
	{
	  Osum[i]/=NAV;
	  if (VERBOSE)
	    printf("i=%i osum=%f\n",i,Osum[i]);
	}
      
      //P5totsum/=(xygrid.gridsize.x[0]/2)*(TAKEN-SKIP);
      
      //printf("volume =%f %i %i \n",xygrid.SIZE*pow(AT,mygrid.d),xygrid.SIZE,mygrid.d);
      
      
      if (VERBOSE)
	{
	  printf("Ut=%f Us=%f vol=%i, Taken=%i skip=%i \n",Utsum,Ussum,vol,TAKEN,SKIP);
	  //for (int x0=0;x0<xygrid.gridsize.x[0]/2;x0++)
	  //  printf("x0=%i P5=%f\n",x0,P5sum[x0]);
	}
      

      //printf("plaqsum=%f, vol=%i, Taken=%i skip=%i\n",Plaqsum,vol,TAKEN,SKIP);
      
      double dW=0;
      double dP=0;
      double dPt=0;
      double dZ=0;
      double dO[largeD];

      
      double dUt=0,dUs=0;
      double dP5sum[xygrid.gridsize.x[0]/2];
      int dWhist[bins],dPhist[bins];
      //double dW2=0;

      
      if (GET_EV_DIST)
	{
	  
	  for (int i=0;i<bins;i++)
	    {
	      dWhist[i]=0;
	      dPhist[i]=0;
	    }
	}
      
      
      for (int x0=0;x0<xygrid.gridsize.x[0]/2;x0++)
	 dP5sum[x0]=0;

      
      for (int conf=SKIP;conf<TAKEN;conf+=STRIDE)//second pass
	{
	  sprintf(buffer,"%s/U-N%i-Nt%i-Ns%i-BETA%.2f-DT%.4f-U%.4f-CT%i-%i.conf",outdir,NC,mygrid.x[0],mygrid.x[1],BETA,dt,UPDATE2*dt,CONFTAKE,conf);
	  loadconf(buffer);
	  if (largeD>0)
	    PolyakovLoop(evs);
	  if (largeD>1)
	    WilsonLoop(evsW,1);
	  double Ut,Us;
	  if (GET_P5_CORR)
	    scalarLoopcorr(P5);
	  
	  PlaqExp(&Ut,&Us);

	  double Otemp[largeD];
	  for (int i=0;i<largeD;i++)
	    Otemp[i]=0;

	  if (largeD>0)
	    Otemp[0]=evs[0];
	  if (largeD>1)
	    Otemp[1]=evsW[0];

	      
	  //double em=EnergyMag();

	  //if (VERBOSE)
	  //  printf("second pass conf=%i Wt=%f Wx=%f Plaq=%f\n",conf,evs[0],evsW[0],em*BETA/vol);
	  dW+=(evsW[0]-Wssum)*(evsW[0]-Wssum);
	  dPt+=(evs[0]-Psum)*(evs[0]-Psum);

	  double temp2=0;
	  if (evs[0]>0.5)
	    temp2+=1.0;
	  if (evsW[0]>0.5)
	    temp2+=1.0;

	  //double temp=evs[0]*evs[0]-W2sum;
	  //temp+=evsW[0]*evsW[0];
	  if (largeD>2)
	    {
	      for (int dir=2;dir<largeD;dir++)
		{
		  WilsonLoop(evsW,dir);
		  //temp+=evsW[0]*evsW[0];
		  //if (evsW[0]>0.5)
		  //  temp2+=1.0;
		  Otemp[dir]=evsW[0];
		}
	    }
	    
	  std::vector<double> myvector (Otemp, Otemp+largeD);
	  std::sort (myvector.begin(), myvector.begin()+largeD);

	  for (int i=0;i<largeD;i++)
	    {
	      //printf("i=%i v=%f\n",i,myvector[i]);
	      dO[i]+=(myvector[i]-Osum[i])*(myvector[i]-Osum[i]);
	    }
	  
	  //dW2+=temp*temp;
	  //dZ+=(temp2-Zsum)*(temp2-Zsum);

	  dUt+=(Ut*BETA/vol-Utsum)*(Ut*BETA/vol-Utsum);
	  dUs+=(Us*BETA/vol-Ussum)*(Us*BETA/vol-Ussum);
	  
	  

	  
	  for (int x0=0;x0<xygrid.gridsize.x[0]/2;x0++)
	    {
	      dP5sum[x0]+=0.5*((P5[x0]-P5sum[x0])*(P5[0]-P5sum[0])+(P5[xygrid.gridsize.x[0]-x0-1]-P5sum[xygrid.gridsize.x[0]-x0-1])*(P5[0]-P5sum[0]));
	      //printf("x0=%i p5=%f s=%f ss=%.8g\n",x0,P5[x0],P5sum[x0],dP5sum[x0]);
	    }

	  printf("Conf %i second pass done!\n",conf);
	  
	}
      dW/=NAV*NAV;
      //dW2/=(TAKEN-SKIP)*(TAKEN-SKIP);
      //dP/=(TAKEN-SKIP)*(TAKEN-SKIP);
      dPt/=NAV*NAV;
      dUt/=NAV*NAV;
      dUs/= NAV*NAV;
      //dZ/=(TAKEN-SKIP)*(TAKEN-SKIP);

      for (int i=0;i<largeD;i++)
	dO[i]/=NAV*NAV;
	
      for (int x0=0;x0<xygrid.gridsize.x[0]/2;x0++)
	dP5sum[x0]/=NAV*NAV;
      
      rav << NAV << "\t";
      rav << Psum << "\t";
      rav << sqrt(dPt) << "\t";
      rav << Wssum << "\t";
      rav << sqrt(dW) << "\t";
      rav << Utsum << "\t";
      rav << sqrt(dUt) << "\t";
      rav << Ussum << "\t";
      rav << sqrt(dUs) << "\t";
      rav << mygrid.x[0] << "\t";
      rav << mygrid.x[1] << "\t";      
      rav << BETA << "\t";
      rav << endl;

      rav.close();

      
      if (GET_ORDERED_WILSON)
	{
	  //W2 << "#NC" << "\t" << "BETA" << "\t" << "Ns" <<"\t"<< "W2" << "\t" << "dW2" << "\t" << "Z \t dZ \t ordered Ws\t\t\t\t"  << "conf" << "\t" << "SKIP" << endl;
	  W2 << "#NC" << "\t" << "BETA" << "\t" << "Ns" <<"\t"<< "ordered Ws\t\t\t\t"  << "tot conf" << "\t" << "SKIP" << endl;
	  //W2 << NC    << "\t" << BETA << "\t" << xygrid.gridsize.x[1] << "\t" << W2sum << "\t" << sqrt(dW2) << "\t" << Zsum << "\t" << sqrt(dZ) << "\t";
	  W2 << NC    << "\t" << BETA << "\t" << xygrid.gridsize.x[1] << "\t";
	

	  for (int i=0;i<largeD;i++)
	    W2 << Osum[i] << "\t" << sqrt(dO[i]) <<"\t";
	  
	  W2 << TAKEN<< "\t" << SKIP << endl;
	}
      
      if (GET_P5_CORR)
	{
	  for (int x0=0;x0<xygrid.gridsize.x[0]/2;x0++)
	    {
	      P5corr << x0 << "\t";
	      P5corr << dP5sum[x0] << "\t";
	      //P5corr << sqrt(dP5sum[x0]) << "\t";
	  P5corr << endl;
	    }
	}

      if (GET_EV_DIST)
	{
	  for (int i=0;i<bins;i++)
	    {
	      double dummydouble=-M_PI+(i+0.5)*width;
	      Hist << dummydouble << "\t" << Phist[i]/(1.0*NAV);
	      Hist << "\t" << Whist[i]/(1.0*NAV);
	      Hist << endl;
	    }
	}

      if (GET_ORDERED_WILSON)
	W2.close();
      if(GET_P5_CORR)
	P5corr.close();
      if (GET_EV_DIST)
	Hist.close();
      
      printf("Nt=%i Ns=%i BETA=%.2f <Wt>=%f+-%f <Wx>=%f+-%f Ut=%f+-%f Us=%f+-%f\n",mygrid.x[0],mygrid.x[1],BETA,Psum,sqrt(dPt),Wssum,sqrt(dW),Utsum,sqrt(dUt),Ussum,sqrt(dUs));
       
      
      free_memory();
      delete [] gens;
    }

if (!ABORT)
     printf("Finished successfully\n");
   else
     printf("Finished with errors\n");
}
