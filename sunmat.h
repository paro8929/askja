#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>


//this can be an algebra matrix or a lie matrix, it does not matter for the 
//implementation. 

class SUNmatrix {

 public:
  int n; //number of rows and columns
  gsl_matrix_complex *M;
  

  SUNmatrix(int myn)
    {
      n=myn;
      M=gsl_matrix_complex_calloc(n,n); //allocate memory and initialize with zeros
    }

  SUNmatrix()
    {      
    }

  ~SUNmatrix()
    {
      //for debugging
      //printf("trying to free %zx \n",M);
      gsl_matrix_complex_free(M);
    }

  void SUNmatrix_init(int nc)
  {
    n=nc;
    M=gsl_matrix_complex_calloc(n,n); //allocate memory and initialize with zeros
  }

  void set(int i,int j,gsl_complex cc)
  {
    gsl_matrix_complex_set(M,i,j,cc);
  }
  
  SUNmatrix& operator*=(const double &d)
  {
    gsl_complex cd;
    GSL_SET_COMPLEX(&cd, d, 0);
    gsl_matrix_complex_scale(M,cd);
    return *this;
  }

  SUNmatrix& operator/=(const double &d)
  {
    gsl_complex cd;
    GSL_SET_COMPLEX(&cd, 1./d, 0);
    gsl_matrix_complex_scale(M,cd);
    return *this;
  }

  SUNmatrix& operator*=(const gsl_complex &d)
  {
    gsl_matrix_complex_scale(M,d);
    return *this;
  }

  SUNmatrix& operator/=(const gsl_complex &d)
  {
    gsl_matrix_complex_scale(M,d);
    return *this;
  }

  SUNmatrix& operator=(const SUNmatrix & rhs)
  {
    //printf("why %i %i %i\n",int(M->tda),int(M->size1),int(M->size2));
    //printf("why2 %i %i %i\n",int(rhs.M->tda),int(rhs.M->size1),int(rhs.M->size2));
    gsl_matrix_complex_memcpy (M,rhs.M); //dont do anything
    return *this;
  }

  SUNmatrix operator*(const SUNmatrix & rhs)
  {
    gsl_complex dummy;
    GSL_SET_COMPLEX(&dummy, 1, 0);
    SUNmatrix temp(n);    

    if(n!=rhs.n)
      printf("Error in Matrix multiplication: %i != %i\n",n,rhs.n);
    
    //printf("checking M: %i %i %i %i\n",n,int(M->size1),int(M->size2),int(M->tda));
    //printf("checking rhs: %i %i %i %i\n",rhs.n,int(rhs.M->size1),int(rhs.M->size2),int(rhs.M->tda));
    //printf("checking temp: %i %i %i %i\n",temp.n,int(temp.M->size1),int(temp.M->size2),int(temp.M->tda));

    //this seems to throw odd errors here and there
    int s=gsl_blas_zgemm (CblasNoTrans,CblasNoTrans,dummy,M,rhs.M,dummy,temp.M);
  
    /*
    for (int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	{
	  GSL_SET_COMPLEX(&dummy, 0, 0);
	  for(int k=0;k<n;k++)
	    dummy=gsl_complex_add(dummy,gsl_complex_mul(gsl_matrix_complex_get(M,i,k),
							gsl_matrix_complex_get(rhs.M,k,j)));	    
	  gsl_matrix_complex_set(temp.M,i,j,dummy);
	}
    //temp.M->size1=M->size1;
	  
    printf("last checking temp: %i %i %i %i\n",temp.n,int(temp.M->size1),int(temp.M->size2),int(temp.M->tda));
    printf("last checking M: %i %i %i %i\n",n,int(M->size1),int(M->size2),int(M->tda));
    printf("last checking rhs: %i %i %i %i\n",rhs.n,int(rhs.M->size1),int(rhs.M->size2),int(rhs.M->tda));*/
	
    return temp;
  }

  SUNmatrix operator+(const SUNmatrix& rhs)
  {
    gsl_complex dummy;
    int enni=rhs.n;
    SUNmatrix temp(enni);    
    for (int i=0;i<enni;i++)
      for(int j=0;j<enni;j++)
	{
	  dummy=gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_matrix_complex_get(rhs.M,i,j));
	  gsl_matrix_complex_set(temp.M,i,j,dummy);
	}
    return temp;
  }

  SUNmatrix operator-(const SUNmatrix& rhs)
  {
    gsl_complex dummy;
    SUNmatrix temp(n);    
    for (int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	{
	  dummy=gsl_complex_sub(gsl_matrix_complex_get(M,i,j),gsl_matrix_complex_get(rhs.M,i,j));
	  gsl_matrix_complex_set(temp.M,i,j,dummy);
	}
    return temp;
  }

  SUNmatrix operator/(const double& d)
  {
    gsl_complex dummy;
    GSL_SET_COMPLEX(&dummy, 1./d, 0);
    SUNmatrix temp(n); 
    gsl_matrix_complex_memcpy (temp.M,M);
    gsl_matrix_complex_scale(temp.M,dummy);
    return temp;
  }

  SUNmatrix operator*(const double& d)
  {
    gsl_complex dummy;
    GSL_SET_COMPLEX(&dummy, d, 0);
    SUNmatrix temp(n); 
    gsl_matrix_complex_memcpy (temp.M,M);
    gsl_matrix_complex_scale(temp.M,dummy);
    return temp;
  }

  SUNmatrix operator*(const gsl_complex & dummy)
  {
    SUNmatrix temp(n); 
    gsl_matrix_complex_memcpy (temp.M,M);
    gsl_matrix_complex_scale(temp.M,dummy);
    return temp;
  }

  

  

  void print()
    {      
      //printf("This is print %i %i %i\n",M->size1,M->size2,M->tda);
      printf("\n");
      for (int i=0;i<n;i++)
	{	 
	  printf("(");
	  for(int j=0;j<n;j++)
	    {
	      //printf("i=%i j=%i\n",i,j);
	      printf("%g+%g i",gsl_matrix_complex_get(M,i,j).dat[0],gsl_matrix_complex_get(M,i,j).dat[1]);
	      if (j!=n-1)
		printf("\t");
	    }
	  printf(")\n");
	}
    }

  

};



//exponent of i times the matrix
SUNmatrix mexpi(const SUNmatrix & abc,int mylevel)
{
  
  gsl_complex dummy;
  GSL_SET_COMPLEX(&dummy, 0, 1);
  SUNmatrix temp(abc.n); 
 
  SUNmatrix rhs(abc.n);
 
  gsl_matrix_complex_set_identity (temp.M);
 
  rhs=abc;
 
  gsl_matrix_complex_scale(rhs.M,dummy);   
 
  switch(mylevel)
    {
    case 0:
      temp=temp;
      break;
    case 1:
      temp=temp+rhs;
      break;
    case 2:
      temp=temp+rhs*(temp+rhs/2.);
      break;
    case 3:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.));
      break;
    case 4:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.)));
      break;
    case 5:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.*(temp+rhs/5.))));
      break;
    case 6:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.*(temp+rhs/5.*(temp+rhs/6.)))));
      break;
    case 7:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.*(temp+rhs/5.*(temp+rhs/6.*(temp+rhs/7.))))));
      break;
    case 8:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.*(temp+rhs/5.*(temp+rhs/6.*(temp+rhs/7.*(temp+rhs/8.)))))));
      break;
    default:
      temp=temp+rhs*(temp+rhs/2.*(temp+rhs/3.*(temp+rhs/4.*(temp+rhs/5.*(temp+rhs/6.*(temp+rhs/7.*(temp+rhs/8.)))))));
    }
  return temp;
}


SUNmatrix Hermitian(const SUNmatrix& abc)
{
  int n=abc.n;
  SUNmatrix temp(n); 
    gsl_complex dummy;
    for (int i=0;i<n;i++)
      for (int j=0;j<n;j++)
	{
	  dummy=gsl_complex_conjugate(gsl_matrix_complex_get(abc.M,i,j));
	  gsl_matrix_complex_set(temp.M,j,i,dummy);
	}

    return temp;
  } 


SUNmatrix transpose(const SUNmatrix& abc)
{
  int n=abc.n;
  SUNmatrix temp(n); 
  
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      gsl_matrix_complex_set(temp.M,j,i,gsl_matrix_complex_get(abc.M,i,j));
  
  
  return temp;
} 

double TrR(const SUNmatrix& abc)
  {
    //printf("funny %i %i %i\n",abc.M->size1,abc.M->size2,abc.M->tda);
    int n=abc.n;
    
     gsl_complex dummy;
     GSL_SET_COMPLEX(&dummy, 0, 0);

     for (int i=0;i<n;i++)
       dummy=gsl_complex_add(dummy,gsl_matrix_complex_get(abc.M,i,i));


     //printf("funny %i %i %i\n",abc.M->size1,abc.M->size2,abc.M->tda);
     return GSL_REAL(dummy);
   }

//imaginary part of trace
double TrI(const SUNmatrix& abc)
  {
    int n=abc.n;
    
     gsl_complex dummy;
     GSL_SET_COMPLEX(&dummy, 0, 0);

     for (int i=0;i<n;i++)
       dummy=gsl_complex_add(dummy,gsl_matrix_complex_get(abc.M,i,i));

     return GSL_IMAG(dummy);
   }

gsl_complex Tr(const SUNmatrix& abc)
{
  //printf("funny %i %i %i\n",abc.M->size1,abc.M->size2,abc.M->tda);
  int n=abc.n;
  
  gsl_complex dummy;
  GSL_SET_COMPLEX(&dummy, 0, 0);
  
  for (int i=0;i<n;i++)
    dummy=gsl_complex_add(dummy,gsl_matrix_complex_get(abc.M,i,i));
  
  
  
  //printf("funny %i %i %i\n",abc.M->size1,abc.M->size2,abc.M->tda);
  return dummy;
}





double tr2(const SUNmatrix& abc)
{
  int n=abc.n;
  SUNmatrix temp(n);
  temp=abc;

  return TrR(temp*temp);
}



//
// transport left computes Udagger.M.U 
//
SUNmatrix transportLeft(SUNmatrix &m,SUNmatrix &u)
{
  int n=m.n;
  SUNmatrix temp(n); 

  temp=Hermitian(u)*m*u;
  return temp;
}

//computes SUN Matrix from Lie Matrix
SUNmatrix Adjoint(const SUNmatrix & abc)
{
  int n=abc.n;
  SUNmatrix temp(n),dummy(n),unity(n);
  gsl_complex minusi;
  GSL_SET_COMPLEX(&minusi, 0, -1);
  

  gsl_matrix_complex_set_identity (unity.M);

  dummy=abc;  
  dummy=dummy-Hermitian(abc);
  temp=(dummy-unity*Tr(dummy)/n)/2.;
  temp=temp*minusi;

  return temp;  
} 
