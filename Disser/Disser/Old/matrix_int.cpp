#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
//на первом месте - возвращаемое значение, дальше - сначала объект, за ним его размерность
//умножение
void multiply(int** C,int const*const*A,const int m,const int n1, int const*const*B,const int n2,const int k)
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

void multiply(int** C,const int*A, const int n1, const int *B, const int n2)
{
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;
		 for (int i=0;i<n;i++)
		 	 for (int j=0;j<n;j++)	  		 
	 	          C[i][j]=A[i]*B[j];
	 	          
		 return ;
}

void multiply(int* C,int const*const*A, const int m,const int n1, const int *B,const int n2)
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

void multiply(int* C,const int*B,const int n2, int const*const*A, const int n1,const int m)
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

void multiply(int* a,const int *v,const int n,const int c)
{ 	 
	   	for(int i=0;i<n;i++)
		    a[i]=v[i]*c;
        
        return ;
}

void multiply(int** C,int const*const*A,const int m,const int n,const int c)
{ 	
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]*c;
		 	 	 
		 return ;
}

void multiply(int*a,const int c,const int *v,const int n)
{ 	 
		for(int i=0;i<n;i++)
		    a[i]=v[i]*c;
        
        return;
}

void multiply(int** C,const int c,int const*const*A,const int m,const int n)
{
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]*c;
		 	 	 
		 return ;
}
//*=
void multiply(int**A, int m,int n1,int const*const*B,int n2,int k)
{
 	 
 		 if (n1!=n2)
 		 	return ;
 		 if (n1!=k)
 		 	return ;	
 		 int n=n1;
 		 int **C=new int*[m];
		 for (int i=0;i<m;i++)
		 {
		  	 C[i]=new int[n];
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



 void multiply(int *v,const int n,const int c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]=v[i]*c;
        
        return ;
}

void multiply(int **A,const int m,const int n,const int c)
{ 	 		 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]=A[i][j]*c;
			 
		 return ;
}


//сложение
void plus(int** C,int const*const*A,const int m1,const int n1, int const*const*B,const int m2,const int n2)
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

void plus(int* C,const int*A, const int n1, const int *B, const int n2)
{
 		 if (n1!=n2)
 		 	return;
 		 	
 		 int n=n1;
		 
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          C[i]=A[i]+B[i];
	 	          
		 return;
}

void plus(int* a,const int *v,const int n,const int c)
{
		for(int i=0;i<n;i++)
		    a[i]=v[i]+c;
        
        return ;
}

void plus(int** C,int const*const*A,const int m,const int n,const int c)
{ 
		 for (int i=0;i<n;i++)
		 	 C[i]=new int[n];
		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]+c;
		 	 	 
		 return ;
}

void plusc(int* a,const int c,const int *v,const int n)
{ 
		for(int i=0;i<n;i++)
		    a[i]=v[i]+c;
        
        return ;
}

void plus(int** C,const int c,int const*const*A,const int m, const int n)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]+c;
		 	 	 
		 return;
}
//+=

void plus( int**A,const int m1,const int n1, int const*const*B,const int m2,const int n2)
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

void  plus( int*A,const int n1, const int *B,const int n2)
{
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;		 		 
			  
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          A[i]+=B[i];
	 	        
		 return;
}

void  plus(int *v,const int n,const int c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]=v[i]+c;           
        return;
}

void plus(int **A,const int m,const int n,const int c)
{ 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]+=c;
		 
		 	 	 
		 return;
}


//минус
void minus(int** C,int const*const*A,const int m1,const int n1, int const*const*B,const int m2,const int n2)
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

void minus(int* C,const int*A, const int n1, const int *B, const int n2)
{
 		 if (n1!=n2)
 		 	return;
 		 	
 		 int n=n1;
		 
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          C[i]=A[i]-B[i];
	 	          
		 return;
}

void minus(int* a,const int *v,const int n,const int c)
{
		for(int i=0;i<n;i++)
		    a[i]=v[i]-c;
        
        return ;
}

void minus(int** C,int const*const*A,const int m,const int n,const int c)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=A[i][j]-c;
		 	 	 
		 return ;
}

void minusc(int* a,const int c,const int *v,const int n)
{
		for(int i=0;i<n;i++)
		    a[i]=-v[i]+c;
        
        return ;
}

void minus(int** C,const int c,int const*const*A,const int m, const int n)
{		 	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[i][j]=-A[i][j]+c;
		 	 	 
		 return;
}

//-=

void minus( int**A,const int m1,const int n1, int const*const*B,const int m2,const int n2)
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

void  minus( int*A,const int n1, const int *B,const int n2)
{
 		 
 		 if (n1!=n2)
 		 	return ;
 		 	
 		 int n=n1;		 		 
			  
		 for (int i=0;i<n;i++)		 	 	  		 
	 	          A[i]-=B[i];
	 	        
		 return;
}

void  minus(int *v,const int n,const int c)
{ 	 
		for(int i=0;i<n;i++)
		    v[i]-=c;           
        return;
}

void minus(int **A,const int m,const int n,const int c)
{ 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 A[i][j]-=c;
		 
		 	 	 
		 return;
}


//----------------------------------------------------------

void transpose_matrix(int **C,int const*const*A,const int m,const int n)
{	 
		 for (int i=0;i<m;i++)
		 	 for (int j=0;j<n;j++)
		 	 	 C[j][i]=A[i][j];
		 	 	 
		 return ;
}

int determinant(int const*const*A,const int m,const int n)
{
	   if (m!=n||m!=2)
	   	  return 0;
 	   
 	   return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}

double norma(const int *v,const int n)
{
		double a=0;
		for(int i=0;i<n;i++)
		   a+=v[i]*v[i];
		return sqrt(a);
}

int scalar_product(const int *v1,const int n1,const int *v2,const int n2)
{
 		 if (n1!=n2)
 		 	return 0;
 		 	
 		 int n=n1;
		 int C=0;
		 for (int i=0;i<n;i++)
		 	 C+=v1[i]*v2[i];	
	 	          
		 return C;
}

void stroka(int*C,int const*const*A,const int m,const int n,const int i)
{
 	   for (int j=0;j<n;j++)
		 	 	 C[j]=A[i][j];
		 	 	 
       return ;
}

void stolbec(int* C,int const*const*A,const int m,const int n,const int i)
{
 	   for (int j=0;j<m;j++)
		 	 	 C[j]=A[j][i];
		 	 	 
       return;
}

void put_stroka(int**A,const int m,const int n1,const int *C,const int n2,const int i)
{
	   if(n1!=n2)
	     return;
	     
       int n=n1;
	   
 	   for (int j=0;j<n;j++)
		 	 	 A[i][j]=C[j];
		 	 	 
       return ;
}

void put_stolbec(int**A,const int n1,const int m,const int *C,const int n2,const int i)
{
	   if(n1!=n2)
	     return;
	     
       int n=n1;
	   
 	   for (int j=0;j<n;j++)
		 	 	 A[j][i]=C[j];
		 	 	 
       return ;
}
//
int** rzeros(int**A,const int m,const int n)
{
 	 A=new int*[m];
 	 for (int i=0;i<m;i++)
 	 {
	  	 A[i]=new int[n];
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=0;
     }
     return A;
}

int* rzeros(int*A,const int n)
{
 	 A=new int[n];
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=0;
     
     return A;
}
//запись нулей в массив и вектор
void zeros(int**A,const int n,const int m)
{
 	
 	 for (int i=0;i<m;i++)
 	 {
	  	
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=0;
     }
     return;
}

void zeros(int*A,const int n)
{
 	 
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=0;
     
     return;
}

//создание массива и вектора и запись в них единиц
int** rones(int**A,const int m,const int n)
{
 	 A=new int*[m];
 	 for (int i=0;i<m;i++)
 	 {
	  	 A[i]=new int[n];
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=1;
     }
     return A;
}

int* rones(int*A,const int n)
{
 	 A=new int[n];
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=1;
     
     return A;
}

//запись единиц в массив и вектор
void ones(int**A,const int m,const int n)
{
 	
 	 for (int i=0;i<m;i++)
 	 {
	  	
 	 	 for (int j=0;j<n;j++)
 	 	 	 A[i][j]=1;
     }
     return;
}

void ones(int*A,const int n)
{
 	 
 	 for (int j=0;j<n;j++)
 	 	 	 A[j]=1;
     
     return;
}
//печать

void print(int const*v,int n)
{
 	 for (int i=0;i<n;i++)
 	 	 cout<<v[i]<<"\n";
 	 	 
 	 return;
}

void print(int const*const*A,int m,int n)
{
 	 for (int i=0;i<m;i++)
 	 {
 	 	 for (int j=0;j<n;j++)
 	 	 	 cout<<A[i][j]<<" ";
 	 	 cout<<"\n";
     }	 
 	 	 return;
}

void print_f(int const*v,int n,const char*s)
{
 	 FILE *f;
 	 f=fopen(s,"w");
 	 fprintf(f,"%d\n",n);
 	 for (int i=0;i<n;i++)
 	 	 fprintf(f,"%d\n",v[i]);
 	 
	  fclose(f);	 
 	 	 return;
}

void print_f(int const*const*A,int m,int n,const char*s)
{
 	 FILE *f;
 	 f=fopen(s,"w");
 	 fprintf(f,"%d\n",m);
 	 fprintf(f,"%d\n",n);
 	 for (int i=0;i<m;i++)
 	 {
 	 	 for (int j=0;j<n;j++)
 	 	 	fprintf(f,"%d ",A[i][j]);
 	 	 fprintf(f,"\n");
     }
 	 fclose(f);	 
 	 	 return;
}
//удаление


void destroy(int **A,int m)
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

void destroy1(int **A,int m)
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
