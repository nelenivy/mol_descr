#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
//на первом месте - возвращаемое значение, дальше - сначала объект, за ним его размерность
//умножение

 	  
void multiply(double** C,double const*const*A,const int m,const int n1, double const*const*B,const int n2,const int k)
{
 		
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<k;j++)
	  		 {
			  	 	  C[i][j]=0;
					  for (int l=0;l<n;l++)
					      C[i][j]+=A[i][l]*B[l][j];
			  }	
			  
		 return;
}

void multiply(double** C,const double*A, const int n1, const double *B, const int n2)
{
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;
		 for (int i=0;i<n;i++)
		 	 for (int j=0;j<n;j++)	  		 
	 	          C[i][j]=A[i]*B[j];
	 	          
		 return ;
}

void multiply(double* C,double const*const*A, const int m,const int n1, const double *B,const int n2)
{
 		 if (n1!=n2)
 		 	return;
 		 	
 		 int n=n1;
		  
		 for (int i=0;i<m;i++)		 	 
		 {
	         C[i]=0;
             for (int l=0;l<n;l++)
					      C[i]+=A[i][l]*B[l];
	     }	
			  
		 return;
}

void multiply(double* C,const double*B,const int n2, double const*const*A, const int n1,const int m)
{
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;
		  
		 for (int i=0;i<m;i++)		 	 
		 {
	         C[i]=0;
             for (int l=0;l<n;l++)
					      C[i]+=A[l][i]*B[l];
	     }	
			  
		 return ;
}

void multiply(double* a,const double *v,const int n,const double c)
{ 	 
	   	for(int i=0;i<n;i++)
		    a[i]=v[i]*c;
        
        return ;
}

void multiply(double** C,double const*const*A,const int m,const int n,const double c)
{ 	
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]*c;
		 	 	 
		 return ;
}

void multiply(double*a,const double c,const double *v,const int n)
{ 	 
		for(int i=0;i<n;i++)
		    a[i]=v[i]*c;
        
        return;
}

void multiply(double** C,const double c,double const*const*A,const int m,const int n)
{
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]*c;
		 	 	 
		 return ;
}
//*=
void multiply(double**A, int m,int n1,double const*const*B,int n2,int k)
{
 	 
 		 if (n1!=n2)
 		 	return ;
 		 if (n1!=k)
 		 	return ;	
 		 int n=n1;
 		 double **C=new double*[m];
		 for (int i=0;i<m;i++)
		 {
		  	 C[i]=new double[n];
		 	 for (int j=0;j<n;j++)
	  		 {
			  	 	  C[i][j]=0;
					  for (int l=0;l<n;l++)
					      C[i][j]+=A[i][l]*B[l][j];
			  }	
	     }	  
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)	  		 
			  	 	  A[i][j]=C[i][j];
			  	 	  
	     for (int i=0;i<m;i++)
	         delete[] C[i];
	         
		 delete C;	  
		 return;
}



 void multiply(double *v,const int n,const double c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]=v[i]*c;
        
        return ;
}

void multiply(double **A,const int m,const int n,const double c)
{ 	 		 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]=A[i][j]*c;
			 
		 return ;
}


//сложение
void plus(double** C,double const*const*A,const int m1,const int n1, double const*const*B,const int m2,const int n2)
{
 		 if (n1!=n2||m1!=m2)
 		 	return ;
 		 	
 		 int n=n1;
 		 int m=m1;
		 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)	  		 
			  	 	  C[i][j]=A[i][j]+B[i][j];					  
			  
		 return ;
}

void plus(double* C,const double*A, const int n1, const double *B, const int n2)
{
 		 if (n1!=n2)
 		 	return;
 		 	
 		 int n=n1;
		 
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          C[i]=A[i]+B[i];
	 	          
		 return;
}

void plus(double* a,const double *v,const int n,const double c)
{
		for(int i=0;i<n;i++)
		    a[i]=v[i]+c;
        
        return ;
}

void plus(double** C,double const*const*A,const int m,const int n,const double c)
{ 
		 for (int i=0;i<n;i++)
		 	 C[i]=new double[n];
		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]+c;
		 	 	 
		 return ;
}

void plus(double* a,const double c,const double *v,const int n)
{ 
		for(int i=0;i<n;i++)
		    a[i]=v[i]+c;
        
        return ;
}

void plus(double** C,const double c,double const*const*A,const int m, const int n)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]+c;
		 	 	 
		 return;
}
//+=

void plus( double**A,const int m1,const int n1, double const*const*B,const int m2,const int n2)
{
 		 if (n1!=n2||m1!=m2)
 		 	return ;
 		 	
 		 int n=n1;
 		 int m=m1;		 
			  
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)	  		 
			  	 	  A[i][j]+=B[i][j];					  
		 
		 return;
}

void  plus( double*A,const int n1, const double *B,const int n2)
{
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;		 		 
			  
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          A[i]+=B[i];
	 	        
		 return;
}

void  plus(double *v,const int n,const double c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]=v[i]+c;           
        return;
}

void plus(double **A,const int m,const int n,const double c)
{ 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]+=c;
		 
		 	 	 
		 return;
}


//минус
void minus(double** C,double const*const*A,const int m1,const int n1, double const*const*B,const int m2,const int n2)
{
 		 if (n1!=n2||m1!=m2)
 		 	return ;
 		 	
 		 int n=n1;
 		 int m=m1;
			  
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)	  		 
			  	 	  C[i][j]=A[i][j]-B[i][j];					  
			  
		 return ;
}

void minus(double* C,const double*A, const int n1, const double *B, const int n2)
{
 		 if (n1!=n2)
 		 	return;
 		 	
 		 int n=n1;
		 
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          C[i]=A[i]-B[i];
	 	          
		 return;
}

void minus(double* a,const double *v,const int n,const double c)
{
		for(int i=0;i<n;i++)
		    a[i]=v[i]-c;
        
        return ;
}

void minus(double** C,double const*const*A,const int m,const int n,const double c)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]-c;
		 	 	 
		 return ;
}

void minus(double* a,const double c,const double *v,const int n)
{
		for(int i=0;i<n;i++)
		    a[i]=-v[i]+c;
        
        return ;
}

void minus(double** C,const double c,double const*const*A,const int m, const int n)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=-A[i][j]+c;
		 	 	 
		 return;
}

//-=

void minus( double**A,const int m1,const int n1, double const*const*B,const int m2,const int n2)
{
 		 if (n1!=n2||m1!=m2)
 		 	return ;
 		 	
 		 int n=n1;
 		 int m=m1;		 
			  
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)	  		 
			  	 	  A[i][j]-=B[i][j];					  
		 
		 return;
}

void  minus( double*A,const int n1, const double *B,const int n2)
{
 		 
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;		 		 
			  
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          A[i]-=B[i];
	 	        
		 return;
}

void  minus(double *v,const int n,const double c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]-=c;           
        return;
}

void minus(double **A,const int m,const int n,const double c)
{ 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]-=c;
		 
		 	 	 
		 return;
}


//----------------------------------------------------------

void transpose_matrix(double **C,double const*const*A,const int m,const int n)
{	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[j][i]=A[i][j];
		 	 	 
		 return ;
}

double determinant(double const*const*A,const int m,const int n)
{
	   if (m!=n||m!=2)
	   	  return 0;
 	   
 	   return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}

double norma(const double *v,const int n)
{
		double a=0;
		for(int i=0;i<n;i++)
		   a+=v[i]*v[i];
		return sqrt(a);
}

double scalar_product(const double *v1,const int n1,const double *v2,const int n2)
{
 		 if (n1!=n2)
 		 	return 0;
 		 	
 		 int n=n1;
		 double C=0;
		 for (int i=0;i<n;i++)
		 	 C+=v1[i]*v2[i];	
	 	          
		 return C;
}

void stroka(double*C,double const*const*A,const int m,const int n,const int i)
{
 	   for (int j=0;j<n;j++)
		 	 	 C[j]=A[i][j];
		 	 	 
       return ;
}

void stolbec(double* C,double const*const*A,const int m,const int n,const int i)
{
 	   for (int j=0;j<m;j++)
		 	 	 C[j]=A[j][i];
		 	 	 
       return;
}

void put_stroka(double**A,const int m,const int n1,const double *C,const int n2,const int i)
{
	   if(n1!=n2)
	     return;
	     
       int n=n1;
	   
 	   for (int j=0;j<n;j++)
		 	 	 A[i][j]=C[j];
		 	 	 
       return ;
}

void put_stolbec(double**A,const int n1,const int m,const double *C,const int n2,const int i)
{
	   if(n1!=n2)
	     return;
	     
       int n=n1;
	   
 	   for (int j=0;j<n;j++)
		 	 	 A[j][i]=C[j];
		 	 	 
       return ;
}
//обратная матрица
double ** obratnaya_matrica(double **E,double const*const*A,const int m,const int n)
{
       if (m!=n)
       	  return 0;
       	  
 	   double **temp_m=new double*[m],*temp_v=new double[m],*temp_v1;
 	   for (int i=0;i<n;i++) 	   
		 	 temp_m[i]=new double[m];
		 	 
		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 temp_m[i][j]=A[i][j];
	
       for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	     if (i!=j)
		 	 	 	E[i][j]=0;
	 	 	     else
	 	 	        E[i][j]=1;
	 	 	        
	 	for (int i=0;i<n;i++)
	 	{
		 	double max=0;
		 	int n_max=0;
		    for (int j=i;j<n;j++) 	        
				if (fabs(temp_m[i][j])>max)
	   			{
				   max=fabs(temp_m[i][j]);
				   n_max=j;
 			    }
 			    
			temp_v1=temp_m[i];
			temp_m[i]=temp_m[n_max];
			temp_m[n_max]=temp_v1;
			temp_v1=E[i];
			E[i]=E[n_max];
			E[n_max]=temp_v1;
			
			multiply(temp_m[i],m,1/temp_m[i][i]);			
			multiply(E[i],m,1/temp_m[i][i]);
			
			for (int j=0;j<n;j++)
			    if (i!=j)
			    {
	   		       stroka(temp_v,temp_m,m,m,i);
	   		       multiply(temp_v,m,temp_m[j][i]);
				   minus(temp_m[j],m,temp_v,m);
			       stroka(temp_v,E,m,m,i);
	   		       multiply(temp_v,m,temp_m[j][i]);
				   minus(temp_m[j],m,temp_v,m);
                }
        }
 	    for (int i=0;i<n;i++) 	   
		 	 delete[] temp_m[i];
		 	 
		delete[] temp_m;
		delete[] temp_v;	 
        return E;
}
				
			
void least_squares(double *ans,double const*const* R,const int m,const int n,double const*z,const int m1)
{
          if(m1!=m) return;
          double **temp_m=new double*[n],**temp_m1=new double*[n];
       	  double *r=new double[n];
       	  double **Rt=new double*[n];
       
       	  for (int j=0;j<n;j++)
       	   	  Rt[j]=new double[m];
       	   	  
 	   	  for (int i=0;i<n;i++)
 	   	  {
		 	 temp_m[i]=new double[n];
		 	 temp_m1[i]=new double[n];
  		  }
		  transpose_matrix(Rt,R,m,n);
      	  multiply(r,Rt,n,m,z,m);
      	  multiply(temp_m,Rt,n,m,R,m,n);
      	  obratnaya_matrica(temp_m1,temp_m,n,n);      	  
      	  multiply(ans,temp_m1,n,n,r,n);
      	  
      	  for (int j=0;j<n;j++)
       	   delete[] Rt[j];
       	   
		   delete[] Rt;
		   
		   for (int i=0;i<n;i++)
		   {
				 	delete[] temp_m[i];
				 	delete[] temp_m1[i];
           }
		 
		   delete[] temp_m;
		   delete[] temp_m1;
		   delete[] r;
}	   
//создание массива и вектора и запись в них нулей

double** rzeros(double**A,const int m,const int n)
{
 	 A=new double*[m];
 	 for (int i=0;i<m;i++)
 	 {
	  	 A[i]=new double[n];
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=0;
     }
     
     return A;
}

double* rzeros(double*A,const int n)
{
 	 A=new double[n];
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=0;
     
     return A;
}
//запись нулей в массив и вектор

void zeros(double**A,const int n,const int m)
{
 	 
 	 for (int i=0;i<m;i++)
 	 {
	  	
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=0;
     }
     return;
}

void zeros(double*A,const int n)
{
 	 
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=0;
     
     return;
}
//создание массива и вектора и запись в них единиц

double** rones(double**A,const int m,const int n)
{
 	 A=new double*[m];
 	 for (int i=0;i<m;i++)
 	 {
	  	 A[i]=new double[n];
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=1;
     }
     return A;
}

double* rones(double*A,const int n)
{
 	 A=new double[n];
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=1;
     
     return A;
}
//запись единиц в массив и вектор

void ones(double**A,const int m,const int n)
{
 	 
 	 for (int i=0;i<m;i++)
 	 {
	  	
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=1;
     }
     return;
}

void ones(double*A,const int n)
{
 	 
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=1;
     
     return;
}

//печать

void print(double const*v,int n)
{
 	 for (int i=0;i<n;i++)
 	 	 cout<<v[i]<<"\n";
 	 	 
 	 	 return;
}

void print(double const*const*A,int m,int n)
{
 	 for (int i=0;i<m;i++)
 	 {
 	 	 for (int j=0;j<n;j++)
 	 	 	 cout<<A[i][j]<<" ";
 	 	 cout<<"\n";
     }
 	 	 
 	 	 return;
}

void print_f(double const*v,int n,const char*s)
{
 	 FILE *f;
 	 f=fopen(s,"w");
 	 fprintf(f,"%d\n",n);
 	 for (int i=0;i<n;i++)
 	 	 fprintf(f,"%lf\n",v[i]);
 	 
	  fclose(f);	 
 	 	 return;
}

void print_f(double const*const*A,int m,int n,const char*s)
{
 	 FILE *f;
 	 f=fopen(s,"w");
 	 fprintf(f,"%d\n",m);
 	 fprintf(f,"%d\n",n);
 	 for (int i=0;i<m;i++)
 	 {
 	 	 for (int j=0;j<n;j++)
 	 	 	fprintf(f,"%lf ",A[i][j]);
 	 	 fprintf(f,"\n");
     }
 	 fclose(f);	 
 	 	 return;
}
//удаление


void destroy(double **A,int m)
{
 	 for (int i=0;i<m;i++)
 	 {
 	 	 delete[] A[i];
 	 	 A[i]=0;
     }
 	
	 delete[] A; 
	 A=0;	 
 	 	 return;
}

void destroy1(double **A,int m)
{
 	 for (int i=0;i<m;i++)
 	 {
 	 	 free(A[i]);
 	 	 A[i]=0;
     }
 	
	 free(A); 
	 A=0;	 
 	 	 return;
}
