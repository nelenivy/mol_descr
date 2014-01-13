#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "matrix_double.cpp"
#include "matrix_int.cpp"

void attributes(double const*const*Vertices,const int vn,int const*const* Triangles,const int fn,int **&neighbours_vv,int **&neighbours_tt)
{
 	 int k=0;
 	 int* neighbours_number,**neighbours_vt;
	 neighbours_number=rzeros(neighbours_number,vn);
  	 neighbours_tt=rzeros(neighbours_tt,fn,3);
	 neighbours_vt=new int*[vn];
	 neighbours_vv=new int*[vn];
	 
	 for (int i=0;i<fn;i++)
	 	 for (int j=0;j<3;j++)					
			neighbours_number[Triangles[i][j]]++;
		 
     
	 for (int i=0;i<vn;i++)
  	 {
		 neighbours_vt[i]= new int[neighbours_number[i]+1];
  	 	 neighbours_vt[i][0]=0;  	 	 
	 }	  	  
 	 delete[] neighbours_number;
 	 for (int i=0;i<fn;i++)
	 	 for (int j=0;j<3;j++)
		  {
		   	neighbours_vt[Triangles[i][j]][0]++; 	  					
			neighbours_vt[Triangles[i][j]][neighbours_vt[Triangles[i][j]][0]]=i;
		  }
	
	 for (int i=0;i<vn;i++)
  	 {
  	 	 neighbours_vv[i]= new int[neighbours_vt[i][0]+1];
  	 	 neighbours_vv[i][0]=neighbours_vt[i][0];
	 }	  	  	  
 	 for (int i=0;i<vn;i++)
 	 {
	  	 int m=0;
	  	 int s=0;
	  	 int e=0;	  	 
 	 	 for (int j=1;j<neighbours_vt[i][0];j++)
 	 	 {
		  	 m=j;
			 s=0;  		  	 		  	 
		  	 while ((s!=2)&&(m!=neighbours_vt[i][0]))
		  	 {
                s=0;
                m++;
				               
                for (int r=0;r<3;r++)
                    for (int l=0;l<3;l++)
                        if (Triangles[neighbours_vt[i][j]][r]==Triangles[neighbours_vt[i][m]][l])
                        {   
                       		s++;
                            if(Triangles[neighbours_vt[i][m]][l]!=i)
                                e=Triangles[neighbours_vt[i][m]][l];
						}
			 }
	         k=neighbours_vt[i][m];
			 neighbours_vt[i][m]=neighbours_vt[i][j+1];
             neighbours_vt[i][j+1]=k;
             neighbours_vv[i][j]=e;
	     }
		   
 	     for (int k=0;k<3;k++)
   	          for (int l=0;l<3;l++)
               	   if (Triangles[neighbours_vt[i][neighbours_vt[i][0]]][k]==Triangles[neighbours_vt[i][1]][l]&&Triangles[neighbours_vt[i][1]][l]!=i)
                      e=Triangles[neighbours_vt[i][1]][l];
                      
         neighbours_vv[i][neighbours_vt[i][0]]=e;
      } 
	  
	   int *temp=new int[fn];
	   
	   for (int i=0;i<fn;i++)
	   	   temp[i]=-1;
	   	   
	   for (int i=0;i<vn;i++)
	   {
	   	   
	   	   int s=0,k=0,l=0;
           for (int j=1;j<=neighbours_vt[i][0];j++)
           {
               k=neighbours_vt[i][j];
               if (j>1)
                  l=neighbours_vt[i][j-1];
        	   else
                  l=neighbours_vt[i][neighbours_vt[i][0]];
            
         	   s=0;
               for (int m=0;m<3;m++)
                   if (neighbours_tt[k][m]==l)
                      s++;
                
           	   if (s==0)
           	   {
                    temp[k]++;
                    temp[l]++;
                    neighbours_tt[k][temp[k]]=l;
                    neighbours_tt[l][temp[l]]=k;
			   }
 		   }
	   }
	   destroy(neighbours_vt,vn);
	   delete[] temp;
	   return;
}



void curvature_cubic(double **H,double const*const*Vertices,const int vn,double const*const*Normals,int const*const*neighbours_vv)
{
       double **E,**temp_m;
       double *temp_v,*r1,*r2,*r3,*ans;
 	   
 	   E=rzeros(E,3,3);
 	   temp_m=rzeros(temp_m,3,3);
 	   temp_v=rzeros(temp_v,3);
 	   r1=rzeros(r1,3);
 	   r2=rzeros(r2,3);
 	   r3=rzeros(r3,3);
 	   ans=rzeros(ans,7);
 	   
 	   int max=0;
 	   
 	   for (int i=0;i<vn;i++)
 	   	   if (neighbours_vv[i][0]>max)
 	   	   	  max=neighbours_vv[i][0];
	   		   	  
       double**neighbours_v,**neighbours_n,**R,*z;
       
       neighbours_v=rzeros(neighbours_v,max,3);
       neighbours_n=rzeros(neighbours_n,max,3);
       R=rzeros(R,3*max,7);
       z=rzeros(z,3*max);  
		 	 
       for (int i=0;i<3;i++)
		 	 for (int j=0;j<3;j++)
		 	     if (i!=j)
		 	 	 	E[i][j]=0;
	 	 	     else
	 	 	        E[i][j]=1;
	 
	  	for (int i=0;i<vn;i++) 
        {  
         	multiply(temp_m,Normals[i],3,Normals[i],3);
			minus(temp_m,3,3,E,3,3);
			multiply(temp_m,3,3,-1);
			temp_v[0]=1;
			temp_v[1]=0; 
			temp_v[2]=0;			        
         	multiply(r1,temp_v,3,temp_m,3,3);
         	multiply(r1,3,1/norma(r1,3));
         	stroka(r3,Normals,vn,3,i);
         	
         	r2[0]=r3[1]*r1[2]-r3[2]*r1[1];
         	r2[1]=-(r3[0]*r1[2]-r3[2]*r1[0]);
         	r2[2]=r3[0]*r1[1]-r3[1]*r1[0];;
         	
            put_stolbec(temp_m,3,3,r1,3,0);
            put_stolbec(temp_m,3,3,r2,3,1);
            put_stolbec(temp_m,3,3,r3,3,2);                 
            
         	for (int j=1;j<=neighbours_vv[i][0];j++)
		    {                  
                     minus(temp_v,Vertices[neighbours_vv[i][j]],3,Vertices[i],3);                     
                     multiply(neighbours_v[j-1],temp_v,3,temp_m,3,3);
                     multiply(neighbours_n[j-1],Normals[neighbours_vv[i][j]],3,temp_m,3,3);
		    }
		    
         for(int j=0;j<neighbours_vv[i][0];j++)
         {
             R[j][0]=pow(neighbours_v[j][0],2)/2;
             R[j][1]=neighbours_v[j][0]*neighbours_v[j][1];
             R[j][2]=pow(neighbours_v[j][1],2)/2;
			 R[j][3]=pow(neighbours_v[j][0],3);
             R[j][4]=pow(neighbours_v[j][0],2)*neighbours_v[j][1];
             R[j][5]=neighbours_v[j][1]*pow(neighbours_v[j][1],2);
             R[j][6]=pow(neighbours_v[j][1],3);
            
             R[neighbours_vv[i][0]+j][0]=neighbours_v[j][0];
             R[neighbours_vv[i][0]+j][1]=neighbours_v[j][1];
             R[neighbours_vv[i][0]+j][2]=0;
             R[neighbours_vv[i][0]+j][3]=3*pow(neighbours_v[j][0],2);
             R[neighbours_vv[i][0]+j][4]=2*neighbours_v[j][0]*neighbours_v[j][1];
             R[neighbours_vv[i][0]+j][5]=pow(neighbours_v[j][1],2);
             R[neighbours_vv[i][0]+j][6]=0;
             
             R[j+2*neighbours_vv[i][0]][0]=0;
             R[j+2*neighbours_vv[i][0]][1]=neighbours_v[j][0];
             R[j+2*neighbours_vv[i][0]][2]=neighbours_v[j][1];
             R[j+2*neighbours_vv[i][0]][3]=0;
             R[j+2*neighbours_vv[i][0]][4]=pow(neighbours_v[j][0],2);
             R[j+2*neighbours_vv[i][0]][5]=2*neighbours_v[j][0]*neighbours_v[j][1];
             R[j+2*neighbours_vv[i][0]][6]=3*neighbours_v[j][1];
             
             z[j]=neighbours_v[j][2];
             
             z[neighbours_vv[i][0]+j]=-neighbours_n[j][0]/neighbours_n[j][2];
             z[neighbours_vv[i][0]*2+j]=-neighbours_n[j][1]/neighbours_n[j][2];
         }
         
      	  least_squares(ans,R,3*neighbours_vv[i][0],7,z,3*neighbours_vv[i][0]);
      	  //print(ans,3);
      	  H[i][0]=4*ans[0]*ans[2]-ans[1]*ans[1];
          H[i][1]=ans[0]+ans[2];
	   }
	   
   for (int j=0;j<max;j++)
   {
	   	  delete[] neighbours_v[j];
	   	  delete[] neighbours_n[j];
   }
   for (int j=0;j<3*max;j++)
   	   delete[] R[j];
   	   
   delete[] neighbours_v;
   delete[] neighbours_n; 
   delete[] R;
   
   for (int i=0;i<3;i++)
 	   {
		 	delete[] E[i];
		 	delete[] temp_m[i];
       }
   
   delete[] ans;      		
   delete[] E;
   delete[] temp_m;
   delete[] temp_v;
   delete[] r1;
   delete[] r2;
   delete[] r3;
   delete[] z;	  
   return;
} 
    
void curvature_parabolloid(double **H,double const*const*Vertices,const int vn,double const*const*Normals,int const*const*neighbours_vv)
{
       double **E=new double*[3],**temp_m=new double*[3];
       double *temp_v=new double[3],*r1=new double[3],*r2=new double[3],*r3=new double[3]; 	   
 	   int max=0;
 	   
 	   for (int i=0;i<vn;i++)
 	   	   if (neighbours_vv[i][0]>max)
 	   	   	  max=neighbours_vv[i][0];
	   		   	  
       double**neighbours=new double*[max],**R=new double*[max],*z=new double[max];
       
	   for (int j=0;j<max;j++)
	   {
	   	   neighbours[j]=new double[3];
	   	   R[j]=new double[3];
	   }
	   
       
 	   for (int i=0;i<3;i++)
 	   {
		 	 E[i]=new double[3];
		 	 temp_m[i]=new double[3];
       }
		 	 
       for (int i=0;i<3;i++)
		 	 for (int j=0;j<3;j++)
		 	     if (i!=j)
		 	 	 	E[i][j]=0;
	 	 	     else
	 	 	        E[i][j]=1;
	 
	  	for (int i=0;i<vn;i++) 
        {  
         	multiply(temp_m,Normals[i],3,Normals[i],3);
			minus(temp_m,3,3,E,3,3);
			multiply(temp_m,3,3,-1);
			temp_v[0]=1;
			temp_v[1]=0; 
			temp_v[2]=0;			        
         	multiply(r1,temp_v,3,temp_m,3,3);
         	multiply(r1,3,1/norma(r1,3));
         	stroka(r3,Normals,vn,3,i);
         	
         	r2[0]=r3[1]*r1[2]-r3[2]*r1[1];
         	r2[1]=-(r3[0]*r1[2]-r3[2]*r1[0]);
         	r2[2]=r3[0]*r1[1]-r3[1]*r1[0];;
         	
            put_stolbec(temp_m,3,3,r1,3,0);
            put_stolbec(temp_m,3,3,r2,3,1);
            put_stolbec(temp_m,3,3,r3,3,2);                 
            
         	for (int j=1;j<=neighbours_vv[i][0];j++)
		    {                  
                     minus(temp_v,Vertices[neighbours_vv[i][j]],3,Vertices[i],3);                     
                     multiply(neighbours[j-1],temp_v,3,temp_m,3,3);
		    }
		    
         for(int j=0;j<neighbours_vv[i][0];j++)
         {
             R[j][0]=neighbours[j][0]*neighbours[j][0];
             R[j][1]=neighbours[j][0]*neighbours[j][1];
             R[j][2]=neighbours[j][1]*neighbours[j][1];             
             z[j]=neighbours[j][2];
         }
         
      	  least_squares(r2,R,neighbours_vv[i][0],3,z,neighbours_vv[i][0]);
      	  
      	  H[i][0]=4*r2[0]*r2[2]-r2[1]*r2[1];
          H[i][1]=r2[0]+r2[2];
	   }
	   
   for (int j=0;j<max;j++)
   {
	   	  delete[] neighbours[j];
	   	  delete[] R[j];
   }
   
   delete[] neighbours; 
   delete[] R;
   
   for (int i=0;i<3;i++)
 	   {
		 	delete[] E[i];
		 	delete[] temp_m[i];
       }
         		
   delete[] E;
   delete[] temp_m;
   delete[] temp_v;
   delete[] r1;
   delete[] r2;
   delete[] r3;
   delete[] z;
   return;	  
} 
