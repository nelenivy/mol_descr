#include <iostream.h>
#include <istream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include<math.h>
#include"singular_points.cpp"
int comp (const void*a,const void*b)
{
 	double r=*(double*)a-*(double*)b;
 	 if (r>0)
 	 	return 1;
 	 else if (r==0)
 	 	  	 return 0;
 	      else
 	         return -1;
}

int main ()
{
 	 ifstream fin;
 	 ofstream fout;
 	 char buf[100],filename[100];
 	 char**opt,**param,*mol_package,*mol_prefix,sng[]=".sng",dst[]=".dst",d[]=".d",md[]=".md";
 	 double**points,*temp,*dist_all,*borders;
	  int*labels,**MD;
	  int size=0,size_all=0,usl=0;
 	 int beg,end,i,mol_number,d_number = 0;
 	 FILE*f;
 	 
 	 //for (d_number = 15; d_number < 16; d_number++)
	   //{
 	 temp=rzeros(temp,3);
 	 opt=new char*[4];
     param=new char*[4];
     
 	 for (int j=0;j<4;j++)
 	 {
	  	 opt[j]=new char[12];
         param[j]=new char[1gik jj,vi;[00];
     }
 	 
 	 strcpy(filename,"options.txt");
 	 strcpy(opt[0],"mol_package");//нужен слэш в конце
 	 opt[1]="mol_prefix";
     opt[2]="mol_number";
     opt[3]="d_number";
	  fin.open(filename);
	  cout<<"File with parameters"<<filename<<"\n";
	  if(fin.fail())
	  {
	   	   cerr<<"invalid filename";
	  	   return 0;
      }
      
      while (!fin.eof())
	  {
	   		fin.getline(buf,sizeof(buf));
	   		
	  		for (i=0;i<4;i++)
	  			if (strstr(buf,opt[i])!=0)
	  			   break;
 	 		if (i<4)//остальные параметры не учитываются, они для построения поверхности
	        {	   
			    for (beg=0;beg<100;beg++)
			    	if (buf[beg]=='"')
			    	   break;
	   	   
	   	        beg++;
	   	        
	   	        for (end=beg;end<100;end++)
			    	if (buf[end]=='"')
			    	   break;
			    	   
	   	        end--;
	   	        
	   	        for (int j=beg;j<=end;j++)
	   	        	param[i][j-beg]=buf[j];
	   	        	
	   			param[i][end-beg+1]='\0';
			}
       }
       
       fin.close();
       mol_package=new char[sizeof(param[0])];
       mol_prefix=new char[sizeof(param[1])];
       strcpy(mol_package,param[0]);
       strcpy(mol_prefix,param[1]);
       mol_number=atoi(param[2]);
	   d_number=atoi(param[3]);
	   
	   //cout<<"особые точки посчитаны? (1 - да,0 - нет)\n";
       //cin>>usl;
       usl = 0;
       
	   	   //temp=rzeros(temp,3);
       if (!usl)
		   for (int j=2;j<=mol_number;j++)
		   {	
		   		cout<<j<<"\n" ;  	   
			  singular_points(points,labels,&size,mol_package,mol_prefix,j);
			  cout<<size<<"1\n";
			  size_all+=(size*(size-1))/2;
			  
		  	  filename[0]='\0';
		  	  cout<<"1\n";
			  sprintf(filename,"%s%s%d%s",mol_package,mol_prefix,j,sng);
			  cout<<"File with singular points"<<filename<<"\n";
			  fout.open(filename);
		  	  
		  	  if(fout.fail())
		  	  {
			   	   cerr<<"invalid filename\n";
			  	   return 0;
	      	   }
	      	   
			   fout<<size<<"\n";
			   
			   for (int k=0;k<size;k++)
			   {
			   	   for (int l=0;l<3;l++)
			   	   	   fout<<points[k][l]<<" ";
			   	   	   
		   		   fout<<labels[k]<<"\n";
			   }
			   
			   fout.close();
			   filename[0]='\0';
			   sprintf(filename,"%s%s%d%s",mol_package,mol_prefix,j,dst);
			   cout<<"File with distances"<<filename<<"\n";
			   fout.open(filename);
		  
		  	   if(fout.fail())
		  	   {
		   	   			   cerr<<"invalid filename\n";
	 				  	   return 0;
	      	   }
	      	   
	      	   fout<<(size*(size-1))/2<<"\n";
	      	   
	      	   for (int k=0;k<size-1;k++)
	      	   	   for (int l=k+1;l<size;l++)
	      	   	   {
				   	   if (labels[l]>labels[k])
		               	  fout<<labels[k]<<" "<<labels[l]<<" ";
        	           else
        	           	   fout<<labels[l]<<" "<<labels[k]<<" ";
        	           	   
				   	   minus(temp,points[l],3,points[k],3);			   	   
	      	   	   	   fout<<norma(temp,3)<<"\n";
	      	   	   	   zeros(temp,3);
	               }
	            
				fout.close();
				destroy(points,size);
				delete[] labels;   
	       }
       else 
	       for (int j=1;j<=mol_number;j++)
		   {
		   	   cout<<j<<"\n";
		   	   
		  	  filename[0]='\0';
			  sprintf(filename,"%s%s%d%s",mol_package,mol_prefix,j,dst);
			  cout<<filename<<"\n";
			  f=fopen(filename,"r");
	      	   fscanf(f,"%d",&size);
	      	   size_all+=size;
	      	   fclose(f);
		   }
	   	   
	   dist_all=rzeros(dist_all,size_all);
	   size_all=0;
	   
	   for (int j=1;j<=mol_number;j++)
		   {
		  	  filename[0]='\0';			  
			   
		  	  filename[0]='\0';
			  sprintf(filename,"%s%s%d%s",mol_package,mol_prefix,j,dst);
			  cout<<filename<<"\n";
			  f=fopen(filename,"r");
 	          fscanf(f,"%d",&size);
	      	   
	      	   int l1,l2;
	      	   for (int k=0;k<size;k++)
	      	   	   fscanf(f,"%d %d %lf",&l1,&l2,&dist_all[size_all+k]);
	      	   	   
  	   		   size_all+=size;
	      	   fclose(f);
		   }
		   
		qsort(dist_all,size_all,sizeof(double),(*comp));
		borders=rzeros(borders,d_number-1);
		
		for (int j=0;j<d_number-1;j++)
		 	 borders[j]=dist_all[(int)floor((j+1)*size_all/d_number)];
		 	 
	    filename[0]='\0';
        sprintf(filename,"%s%s%s%d",mol_package,mol_prefix,d,d_number);
        cout<<filename<<"\n";
        fout.open(filename);
		  	  
  	    if (fout.fail())
  	    {
		   	   	 cerr<<"invalid filename";
		  	   	 return 0;
   	    }
	    
	    for (int j=0;j<d_number-1;j++)
 	         fout<<borders[j]<<"\n";
 	    
	    delete[] dist_all;
		MD=rzeros(MD,mol_number,36*d_number);     
        //cin>>i;
        
        for (int j=1;j<=mol_number;j++)
        {
		 	filename[0]='\0';
		    sprintf(filename,"%s%s%d%s",mol_package,mol_prefix,j,dst);
		    cout<<filename<<"\n";
		    f=fopen(filename,"r");
            fscanf(f,"%d",&size);
	      	   
	      	   int l1,l2,i=0;
	      	   double r;
	      	   
			for (int k=0;k<size;k++)
			{
			 	fscanf(f,"%d %d %lf",&l1,&l2,&r);
			 	i=0;
	      	    while ((i<d_number-1)&&(r>borders[i])) i++;	      	    
	      	    MD[j-1][l1-1+6*(l2-1)+36*i]++;
			}
			
			fclose(f);
		}
	filename[0]='\0';
    sprintf(filename,"%s%s%s%d",mol_package,mol_prefix,md,d_number);
    cout<<filename<<"\n";	
	print_f(MD,mol_number,36*d_number,filename);
	
	
	delete[]temp;
	delete[]borders;
	destroy(MD,mol_number);

 	delete[] mol_package;
	delete[] mol_prefix;
	cin>>i;

	return 1;
}
