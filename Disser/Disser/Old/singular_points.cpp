#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "curvature.cpp"
#include "segmentation.cpp"
#include "marking.cpp"
#include "surface_read.cpp"
            

   	        
void singular_points(double**&points,int*&labels,int *size,const char*mol_package,const char*mol_prefix,const int mol_number)
{
 	double **Vertices,**Normals;
    double **H; 	
	int	**Triangles;
	int **neighbours_vv;
	int **neighbours_tt;
	int *colour_matrix,*Segments,*numbers,*centers;
	int max_size=250,i;
	char filename[200],filename1[200],srf[]="srf",ch[]="ch",wrl[]=".wrl",txt[]=".txt";	
	sprintf(filename,"%s%d%s",mol_package,mol_number,wrl);
	cout<<"File with surface"<<filename<<"\n";
	
	surface_read(filename,Vertices, Triangles, Normals);
	cout<<"1\n";
	int vn=_msize(Vertices)/sizeof(Vertices[0]);
	int fn=_msize(Triangles)/sizeof(Triangles[0]);

	H=rzeros(H,vn,2);
	colour_matrix=rzeros(colour_matrix,vn);
	Segments=rzeros(Segments,vn);
	numbers=rzeros(numbers,3);	

	attributes(Vertices,vn,Triangles,fn,neighbours_vv,neighbours_tt);
	
	curvature_cubic(H,Vertices,vn,Normals,neighbours_vv);
	
	for (int i=0;i<vn;i++)
		if (H[i][1]>=0)//средняя кривизна
		   colour_matrix[i]=1;		
	
	curvature_parabolloid(H,Vertices,vn,Normals,neighbours_vv);
	
	
	
	for (int i=0;i<vn;i++)
		if (colour_matrix[i]==0)
		   if (H[i][0]>=0)
		   	  colour_matrix[i]=-1;
 		   else
 		   	   colour_matrix[i]=0;
 		   	   
     
     cout<<"1\n";
 	//1 - впадина,-1 - выпуклый, 0 - седло	
	 segmentation(Segments,numbers,Vertices,vn,Triangles,fn,neighbours_vv,neighbours_tt,colour_matrix,max_size) ;  	   

	 /*print_f(numbers,3,"numbers.txt");
	 print_f(Segments,vn,"Segments.txt");
	 print_f(Vertices,vn,3,"Vertices.txt");
	 print_f(Triangles,fn,3,"Triangles.txt");*/
	 
	 labels=rzeros(labels,numbers[2]);
	 	
	 centers=rzeros(centers,numbers[2]);
	 
	 sprintf(filename1,"%s%d%s",mol_package,mol_number,txt);
	 cout<<"1\n";
	 cout<<"File with charges"<<filename1<<"\n";
	 marking(labels, centers,filename1,Segments,numbers,Vertices,vn);
	 points=rzeros(points,numbers[2],3);
	 
	 for (int i=0;i<numbers[2];i++)
	 	 for (int j=0;j<3;j++)
	 	 	 points[i][j]=Vertices[centers[i]][j];
	
	 *size=numbers[2];	 	 
	destroy1(Vertices,vn); 
	destroy1(Triangles,fn);
	destroy1(Normals,vn);
	destroy(H,vn); 
	
	destroy(neighbours_vv,vn);
	destroy(neighbours_tt,fn);
	delete[]colour_matrix;
	delete[]Segments;
	delete[]numbers;
	delete[]centers;

	//cin>>i;
	return ;
}
