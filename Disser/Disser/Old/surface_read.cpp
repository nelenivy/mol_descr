 #include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include<math.h>

void surface_read(const char*filename, double**&Vertices, int**&Triangles, double**&Normals)
{
 	 Vertices=(double**)malloc(sizeof(double*)) ;
 	 Triangles=(int**)malloc(sizeof(int*));
 	 
 	 int n=1;
 	 char buf[60];
	  cout<<"1\n"; 	 
 	 ifstream f;
 	 
	 f.open(filename);
	  
 	 if(f.fail())
	  {
	   	   cerr<<"invalid filename\n";
	  	   return;
      }
      f.getline(buf,sizeof(buf));
      
      while (strstr(buf, "Triangle")==0&&(!(f.eof())))      
	   		f>>buf;
	   		
      
	  
      while (strstr(buf, "point")==0&&!f.eof())
	   		f.getline(buf,sizeof(buf));
	   		
      int i=0;      
    
	  while(f.rdstate()==0)
	  {
	   	    Vertices[i]=(double*)malloc(3*sizeof(double)) ;
	   		
            for	(int j=0;j<3;j++)					   
	  		    f>>(Vertices[i][j]);
	  		    
		    
	  		f>>buf;    
			i++;
			//cout<<i<<"\n";				
			if (i+1>n)
			{
	           n+=1000;
			   Vertices=(double**)realloc(Vertices,n*sizeof(double*));
            }
      }
       
      int vn=i-1;      
	  Vertices=(double**)realloc(Vertices,vn*sizeof(double*));
	  f.clear();
	  	
      while (strstr(buf, "# Normal definition")==0&&!f.eof())
	   		f.getline(buf,sizeof(buf));
	   		
      i=0;      
      f.getline(buf,sizeof(buf));
      f.getline(buf,sizeof(buf));
      Normals=(double**)malloc(vn*sizeof(double*));
      
	  for (int i=0;i<vn;i++)
	  {
	   		Normals[i]=(double*)malloc(3*sizeof(double));
	   		
            for	(int j=0;j<3;j++)					   
	  		    f>>(Normals[i][j]);
	  		f>>buf; 
      } 
	     
	  f.clear();
	    
	  while (strstr(buf, "IndexedFaceSet")==0&&!f.eof())
	   		f.getline(buf,sizeof(buf));
		   		
      while (strstr(buf, "coordIndex [")==0&&!f.eof())
	   		f.getline(buf,sizeof(buf));
	   		
      i=0;      
      n=1;
	  while(f.rdstate()==0)
	  {
	   					   
	   		Triangles[i]=(int*)malloc(3*sizeof(int)) ;;
	   		char c;
            for	(int j=0;j<3;j++)
			{					   
	  		    f>>(Triangles[i][j]);
	  		    f>>c;
			}
			
	  		f>>buf;    
			i++;			
			if (i+1>n)
			{
	           n+=1000;
			   Triangles=(int**)realloc(Triangles,n*sizeof(int*));
            }
      }
      
      f.clear();
      int fn=i-1;      
	  Triangles=(int**)realloc(Triangles,fn*sizeof(int*));
	  f.close();
	  return;
}	

