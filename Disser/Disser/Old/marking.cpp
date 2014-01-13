#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
void marking(int*labels, int*centers,char*filename,const int* Segments,const int*numbers,double const*const*Vertices,const int vn)
{
 	 double*p,**coord,*charges,*temp;
     int l,*num,**S,k=numbers[2];
     ifstream f;
     cout<<k<<"\n";
     p=rzeros(p,k);
     temp=rzeros(temp,3);
     num=rzeros(num,k);
     S=new int*[k];
     
 	 for (int i=0;i<vn;i++)
 	 	 if (Segments[i]>0)
		     num[Segments[i]-1]++;
	
	for (int i=0;i<k;i++)
		S[i]=rzeros(S[i],num[i]);
		
	zeros(num,k);	
		
	 for (int i=0;i<vn;i++)
 	 	 if (Segments[i]>0)
 	 	 {
		     num[Segments[i]-1]++;	
		     S[Segments[i]-1][num[Segments[i]-1]-1]=i;
	     }
	     
	for (int i=0;i<k;i++)     
		centers[i]=center(Vertices,vn,S[i],num[i]);	  
		
	 delete[] num; 
	 destroy(S,k);
	  
	 f.open(filename);
	 
 	 if(f.fail())
	  {
	   	   cerr<<"invalid filename\n";
	  	   return;
      }
      f>>l;
      
      coord=rzeros(coord,l,3);
      charges=rzeros(charges,l);
      
      for (int i=0;i<l;i++)
      {
	   	  for (int j=0;j<3;j++)
	   	  	  f>>coord[i][j];
	   	  	  
 	  	  f>>charges[i];
	  }
	  
	  f.close();
	  
	  for (int i=0;i<k;i++)
	  	  for (int j=0;j<l;j++)
	  	  {
		   	  //cout<<i-k<<"\n";
	   	      minus(temp,coord[j],3,Vertices[centers[i]],3);
	   	      p[i]+=charges[j]/norma(temp,3);
	   	      zeros(temp,3);
	      }
	      
      //cout<<"1\n";
	  for (int i=0;i<k;i++)
	  {
	   	  
  	  	 	 
  	  	 int j=0;	 
		 while(i>numbers[j])
			j++;
		
		labels[i]=2*j+1;
		if (p[i]>0)
	   	  	 labels[i]+=1;
  	  	 
	  }
	  
	  delete[] p;
	  destroy(coord,l);
	  delete[] charges;
	  delete[]temp;
	  return;
}
