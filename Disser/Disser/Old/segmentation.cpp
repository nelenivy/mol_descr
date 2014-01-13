#include <iostream>
#include <istream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
//#include "matrix_double.cpp"
//#include "matrix_int.cpp"
void grow_region(int *Segments,int *number,const int i,const int k, int const*const* neighbours_tt,int const* list)
{
 	 Segments[i]=k;
 	 (*number)++;
 	 for(int j=0;j<3;j++)
 	 {
	  		 int m=neighbours_tt[i][j];
 	 		 if (list[m]==1&&Segments[m]==0)
 	 		 	grow_region(Segments,number,m,k,neighbours_tt,list);
     }
     return;
}


int segmentation_c(int *Segments,int const*const* neighbours_tt,int const*list,const int fn)//выделяет связные сегменты вершин
{
 	int k=0,*sn;
 	int *exc,l=0;
 	sn=rzeros(sn,fn);
 	for (int i=0;i<fn;i++)
	 	if (list[i]==1&&Segments[i]==0)
	 	{
   		   k++;
	 	   grow_region(Segments,&sn[k-1],i,k,neighbours_tt,list);
        }
   
   
   exc=rzeros(exc,k+1);
   //удаление маленьких сегментов
   for (int i=1;i<=k;i++)
   	   if (sn[i-1]>5)
   	   {
  		   l++;
  		   exc[i]=l;
       }
       
   k=l;
   
   for (int i=0;i<fn;i++)
   	   Segments[i]=exc[Segments[i]];
   	
	delete[] exc;
	delete[] sn;
	//sn=(int*)realloc(sn,k);
	return k;
}
//------------------------------------------------------------------------------   
 //к-средних
int center(double const*const*Vertices,const int vn,int const *list,const int l)
{
	//ищем элемент с минимальной суммой расстояний до среднего
	//арифметического координат элементов кластера
	
	//c-индекс найденного центра 
	//Vertice_Matrix - список координат всех вершин триангуляции
	//list - массив с номерами вершин, входящих в кластер
	double *d,*x,m=0;
	int v_min=0;
	d=rzeros(d,l);
	//ищем элемент с минимальной суммой расстояний до среднего
	//арифметического координат элементов кластеров
	x=rzeros(x,3);
	for (int j=0;j<l;j++)
	     for (int k=0;k<3;k++)
	         x[k]+=Vertices[list[j]][k]/l;
	
	for (int j=0;j<l;j++)
	{
	    for (int k=0;k<3;k++)
	        d[j]+=pow(Vertices[list[j]][k]-x[k],2);
	    d[j]=sqrt(d[j]);
	}
	for (int j=0;j<l;j++)
		if (d[j]>m)
		{
		   m=d[j];
		   v_min=j;
        }
     delete[] x;
     delete[] d;
	return list[v_min];
}
	
 void k_means(int* Segments,const int l,double const*const*Vertices,double const*const*Distances,const int vn,int const*list,const int vn_1)
 {
  	 
	int s=0;
	int *seed,v_min=0;//центры кластера
	int **segments,*sn;
	int r=0;
	int count=0,s_min;
	double d_min;
	
	s=(int)floor(vn_1/l);
	//cout<<"l="<<l<<"\n";
	//среднее количество элементов в кластере
	sn=rzeros(sn,l);
	segments=rzeros(segments,l,vn_1);
	seed=rzeros(seed,l);
	
	for (int k=0;k<l;k++)
	{
	    int j=0;
	    while(j<s)
	    {
	        if(list[r]>0)
	            j++;
	        r++;
	    }
	    seed[k]=r-1;
	}
	
	//количество изменений центров кластеров во время итерации
	s=1;
	//список кластеров с вершинами им принадлежащими
	//в iой сроке на 1 месте количество элементов в кластере, дальше список
	//вершин, принадлежащих кластеру
	
	while(s>0&&count<1000)
	{
	 					  
        zeros(sn,l);
	    count++;
	    s=0;
	    //относим iую вершину к тому кластеру, к центру которого она ближе
	    
	    
	    
	    zeros(sn,l);
	    //zeros(segments,l,vn_1);
	    for (int i=0;i<vn;i++)
	    {
	        if(list[i]>0)
	        {	 
	            d_min=10000;
	            s_min=0;
	            for (int j=0;j<l;j++)
	                if(Distances[i][seed[j]]<d_min)
	                {
	                    d_min=Distances[i][seed[j]];
	                    s_min=j;
					}
	                       
	            sn[s_min]++;
	            segments[s_min][sn[s_min]-1]=i;
            }  
        }
	    //пересчитываем центры вершин
	    for (int i=0;i<l;i++)
	    {
	        v_min=center(Vertices,vn,segments[i],sn[i]);
	            //смотрим, изменился ли центр с предыдущей итерации
	        if(v_min!=seed[i])
	        {
	                s++;
	                seed[i]=v_min;
	        }      
	    }
	}
	
	 for (int i=0;i<l;i++)     
	     for (int j=0;j<sn[i];j++)
	         Segments[segments[i][j]]=i+1;
	 
	 destroy(segments,l);
	 delete[] sn;
	 delete[] seed;
	 return;
}

int k_means_1(int* Segments, int l,double const*const*Vertices,double const*const*Distances,const int vn)
{
//k_means_1 находит компоненты связности среди вершин на
//основе Distance_Matrix

//Segment_Matrix - массив, на iом месте номер кластера, к которому
//принадлежит iая вершина
//l-количество кластеров
//Distance_Matrix-cодержит кратчаюшие расстояния между вершинами, для
//которых list больше 0, если вершины i,j в разных компонентах связанности,
//то Distance_Matrix(i,j)=0

	int count=0;
	int *list_1,*list;
	int k=0;
	int m=0;
	int *S;
	
	S=rzeros(S,vn);
	list_1=rzeros(list_1,vn);
	list=rones(list,vn);
	
	for (int i=0;i<vn;i+f`sz+)
	//для каждой компоненты связности запускаем к-средних
	    if(list[i]>0)
	    {
	        zeros(list_1,vn);
	        zeros(S,vn);
	        list_1[i]=1;
	        list[i]=0;
	        //сount - общее количество кластеров на предыдущем шаге
	        //m - кол-во кластеров на предыдущем шаге
	        count+=m;
	        k=1;
	        //%list_1(j)=1, если jая вершина находится в одной компоненте
	        //%связности с iой
	        //%те вершины, которые мы кластеризуем, убираем из списка list
	        for (int j=0;j<vn;j++)
	            if(Distances[i][j]>0)
	            {
	                list_1[j]=1;
	                list[j]=0;
	                k++;
	            }
    		   //cout<<"k="<<k<<"\n";
	        m=(int)floor(l*k/vn);
	        
	        //%m-число кластеров на этом шаге
	        if(m!=0)
			{    
	            k_means(S,m,Vertices,Distances,vn,list_1,k);
	            //%пересчитываем номера кластеров
	            for (int j=0;j<vn;j++)
	                if(S[j]>0)
	                    S[j]+=count;
	                    
	            plus(Segments,vn,S,vn);
	        }
	    }
	    //cout<<"---------------\n";
	l=count+m;
	
	delete[] S;
	delete[] list_1;
	delete[] list;
	return l;
}


//------------------------------------------------
void min_heapify(double**Q,int*num,const int i,const int l)
{
//%min_heapify - опускает элемент с iым номером вниз по пирамиде Q,если он
//%нарушает упорядоченность пирамиды
//%соответствущим образом изменяет num

//%Q - массив, на 1 месте расстояние до iой вершины, на 2 - номер
//%вершины
//%num(i) - под каким номером iая вершина находится в Q
//%i - номер элемента, который опускать
//%Очередь Q организована в виде пирамиды, l - длина очереди
   int m=-1;
   
   if((2*i<l)&&(Q[2*i][0]<Q[i][0]))
       m=2*i;
   else
       m=i;
   
   if(2*i+1<l&&Q[2*i+1][0]<Q[m][0])
       m=2*i+1;   
   
   if(i!=m)
   {
       num[(int)Q[m][1]]=i;
       num[(int)Q[i][1]]=m;
       Q[l][0]=Q[m][0];
       Q[l][1]=Q[m][1];
       Q[m][0]=Q[i][0];
       Q[m][1]=Q[i][1];
       Q[i][0]=Q[l][0];
       Q[i][1]=Q[l][1];
       min_heapify(Q,num,m,l);      
   }
   return;
}

void extract_min(double*m,double**Q,int *num,int *l)
{
//%extract_min - достать минимальный элемент из очереди, соответствующим
//%образом изменить num и Q

//%m - найденное минимальное значение
//%Q - массив, на 1 месте расстояние до iой вершины, на 2 - номер
//%вершины
//%num(i) - под каким номером iая вершина находится в Q
//%Очередь Q организована в виде пирамиды, l - длина очереди
	m[0]=Q[0][0];
	m[1]=Q[0][1];
	num[(int)Q[0][1]]=-1;
	num[(int)Q[*(l)-1][1]]=0;
	Q[0][0]=Q[*(l)-1][0];
	Q[0][1]=Q[*(l)-1][1];
	(*l)--;
	min_heapify(Q,num,0,*l);
	return;
}

void decrease_key(double**Q,int*num, int i,const double k,const int l)
{
//%decrease_key - присваивает Q(i,1) значение k, если k меньше,
//%соответствующим образом изменяет num и Q
//
//%Q - массив, на 1 месте расстояние до iой вершины, на 2 - номер
//%вершины
//%num(i) - под каким номером iая вершина находится в Q
//%Очередь Q организована в виде пирамиды, len - длина очереди
//%i - номер элемента, ключ которого надо изменить
//%k - новое значение ключа
	
	if(Q[i][0]>k)
	{
	   Q[i][0]=k;
	   int m=(int)floor(i/2);
	   while(m>=0&&Q[i][0]<Q[m][0])
	   {
	       num[(int)Q[i][1]]=m;
	       num[(int)Q[m][1]]=i;
	       Q[l][0]=Q[m][0];
	       Q[l][1]=Q[m][1];
	       Q[m][0]=Q[i][0];
	       Q[m][1]=Q[i][1];
	       Q[i][0]=Q[l][0];
	       Q[i][1]=Q[l][1];
	       i=m;
	       m=(int)floor(i/2);
	   }
    }
	
	return;  
}


void distance_piram(double**Distances,double const*const*Vertices,int const*const*neighbour_matrix_vv,const int vn)
{
//%distance - вычисляет матрицу расстояний 

//%Distances-cодержит кратчаюшие расстояния между вершинами, если вершины i,j в разных компонентах связанности,
//%то Distances(i,j)=0
//%Vertices - список координат вершин
//%neighbour_matrix_vv-iая строка содержит на
//%1ом месте количество элементов в строке, дальше список вершин, соседних к
//%iой, идущие друг за другом
//%list-список вершин, для которых рассчитываем

	int k=0;
	double **adjacency_matrix;
	int *num,len;
	double **Q,*d_min=new double[2];
	int d,l,r;
	
	adjacency_matrix=rzeros( adjacency_matrix,vn,vn);
	num=rzeros(num,vn);
	//
	//zeros(Distances,vn,vn);
	//%вычисляем матрицу смежности 
	
	for (int i=0;i<vn;i++)
	        for (int j=1;j<=neighbour_matrix_vv[i][0];j++)
	        {
	            k=neighbour_matrix_vv[i][j];
	            
	            if(adjacency_matrix[i][k]<=0)
                {
			  		    adjacency_matrix[i][k]=sqrt(pow(Vertices[i][0]-Vertices[k][0],2)+pow(Vertices[i][1]-Vertices[k][1],2)+pow(Vertices[i][2]-Vertices[k][2],2));
	                    adjacency_matrix[k][i]=adjacency_matrix[i][k];                   
                }            
	        }

	//%вычисляем кратчайший путь
	//%Q - очередь вершин
	
	Q=rzeros(Q,vn+1,2);
	
	for (int i=0;i<vn;i++)	    
 	{	 
	        //%если вершина в списке, формируем очередь из вершин,для которых
	        //%list>0, ставим i  на 1 место
	        //%Q - массив, на 1 месте расстояние до iой вершины, на 2 - номер
	        //%вершины
	        //%num(i) - под каким номером iая вершина находится в Q
	        //%Очередь Q организована в виде пирамиды, len - длина очереди
	        
	        for (int j=0;j<vn;j++)	            
            {
	                Q[j][0]=10000;                
	                Q[j][1]=j;
	                num[j]=j;
 			}
			Q[i][1]=0;
            num[0]=i; 
            Q[0][1]=i;
            num[i]=0;
            Q[0][0]=0; 
	        	        
	        len=vn;
	        while(len>0)
			{          
	        //%находим вершину с минимальным расстоянием до iой вершины из
	        //%очереди
	            extract_min(d_min,Q,num,&len);
				d=(int)d_min[1] ; 
				
	        //%пересчитываем расстояния для соседних вершин
	            for (int j=1;j<=neighbour_matrix_vv[d][0];j++)
	            {
	                l=neighbour_matrix_vv[d][j];
	                //cout<<"len="<<len<<" "<<(int)floor(len/2)<<"\n";
	                if(num[l]!=-1)
	                    decrease_key(Q,num,num[l],d_min[0]+adjacency_matrix[l][d],len);
	                    
						
	            } 
	            
	            if(d_min[0]==10000)
	              d_min[0]=0;
	              
	            Distances[i][d]=d_min[0];
	            Distances[d][i]=d_min[0];
	        }
	    }
	    
		delete[] d_min;
		delete[] num;
		
		destroy(Q,vn);
		destroy(adjacency_matrix,vn);
		return;
		
}
//---------------------------------------------

void mapping(int*map,int*map_1,double**Vertices_1,const int vn_1,double const*const*Vertices, int const*list,const int vn)
{
//%перенумеровывает несегментированные вершины и записывает соответствие в
//%map
     int count=-1;
     for (int i=0;i<vn;i++)
	    if(list[i]>0)
	    {
		 			 
	        count++;
	        map[count]=i;
	        map_1[i]=count;
	        
	        for (int j=0;j<3;j++)
	        {
	         	Vertices_1[count][j]=Vertices[i][j];
	        }
	    }
	    
	return;
}

int ** neighbours(int const*const*neighbours_vv,const int vn,int const*map_1,const int vn_1)
{
 	int*numbers,k; 
	numbers=rzeros(numbers,vn_1);
	
	for (int i=0;i<vn;i++)
	    if(map_1[i]>=0)
	        for (int j=1;j<=neighbours_vv[i][0];j++)
	        {
	            k=neighbours_vv[i][j];
	            if(map_1[k]>=0)
	                numbers[map_1[i]]++;
	        }
	
	int** neighbours_vv_1=new int*[vn_1];
	
	for (int i=0;i<vn_1;i++)
	{
	    neighbours_vv_1[i]=new int[numbers[i]+1];
	    neighbours_vv_1[i][0]=0;
	}
	    
	delete[] numbers;
	
	for (int i=0;i<vn;i++)
	    if(map_1[i]>=0)
	        for (int j=1;j<=neighbours_vv[i][0];j++)
	        {
	            k=neighbours_vv[i][j];
	            if(map_1[k]>=0)
	            {
	                neighbours_vv_1[map_1[i]][0]++;
	                neighbours_vv_1[map_1[i]][neighbours_vv_1[map_1[i]][0]]=map_1[k];
	            }
	        }
	        
     return neighbours_vv_1;
 }
	
int segment_cut(int*Segments_v,double const*const*Vertices,const int vn,int const*const*neighbours_vv,const int max_size,const int i,int k1)
{	
			//cout<<"dghdftgh\n";
			double**Distances_1,**Vertices_1;
            int vn_1=0,*Segments,*Segments_1,*S,*map,*map_1,**neighbours_vv_1;
			S=rzeros(S,vn);
			Segments_1=rzeros(Segments_1,vn);
            for (int j=0;j<vn;j++)
                if(Segments_v[j]==i)
                {
                    S[j]=1;
                    vn_1++;
			    }
            Distances_1=rzeros(Distances_1,vn_1,vn_1);  
			Segments=rzeros(Segments,vn_1);    
			map=rzeros(map,vn_1);
     		Vertices_1=rzeros(Vertices_1,vn_1,3);
     		map_1=rzeros(map_1,vn);
     		minus(map_1,vn,-1);
            mapping(map,map_1,Vertices_1,vn_1,Vertices, S,vn);           
            neighbours_vv_1=neighbours(neighbours_vv,vn,map_1,vn_1);
            //%перенумеровывает несегментированные вершины и записывает соответствие в
			distance_piram(Distances_1,Vertices_1,neighbours_vv_1,vn_1);
            int k=(int)(k1/max_size)+1;
            k=k_means_1(Segments, k,Vertices_1,Distances_1,vn_1);
            //%кластеров может получиться меньше, чем к, потому что слишком маленькие компоненты связности не сегментируются
            
            for (int j=0;j<vn_1;j++)
                Segments_1[map[j]]=Segments[j];
            
            for (int j=0;j<vn;j++)
                if(Segments_1[j]>1)
                   Segments_v[j]=k1+Segments_1[j]-1;
              
            k1+=k-1;
            delete[]Segments;
			delete[]Segments_1;
			delete[]S;
			delete[]map;
			delete[]map_1;
			
			for (int j=0;j<vn_1;j++)
			{
			 	delete[] neighbours_vv_1[j];
			 	delete[]Distances_1[j];
                delete[]Vertices_1[j];
			}
			
			delete[] neighbours_vv_1;
	 	    delete[]Distances_1;
            delete[]Vertices_1;
            return k1;
}
//------------------------------
//segmentacia
void segmentation(int*Segments_v,int *numbers,double const*const*Vertices,const int vn,int const*const* Triangles,const int fn,int const*const*neighbours_vv,int const*const*neighbours_tt,int*colour_matrix,const int max_size)
{
 	int k1,k2,k,*list_t,*sv,*Segments_t,*Segments_1;
 	double**Distances;
 	double**Vertices_1;
  	int vn_1=0,*list_v,*map,*map_1,**neighbours_vv_1;
 	Segments_t=rzeros(Segments_t,fn);
	list_t=rzeros(list_t,fn);
	list_v=rzeros(list_v,vn);
	
 	 for (int i=0;i<vn;i++)
 	{
	 		 int s=0;
			 for (int j=1;j<=neighbours_vv[i][0];j++)
			 {
			  	 if (colour_matrix[i]!=colour_matrix[neighbours_vv[i][j]])
				   s++;
             }
			 if (s==neighbours_vv[i][0])
			 	colour_matrix[i]=colour_matrix[neighbours_vv[i][1]];
	}
	//определение типа треугольника
	for (int i=0;i<fn;i++)
	 	if (colour_matrix[Triangles[i][0]]==colour_matrix[Triangles[i][1]]&&colour_matrix[Triangles[i][0]]==colour_matrix[Triangles[i][2]])
		   	list_t[i]=colour_matrix[Triangles[i][0]];
		   
	//сегментация
	k1=segmentation_c(Segments_t,neighbours_tt,list_t,fn);
	
	sv=rzeros(sv,k1);
	//print(Segments_t,fn);
	for (int i=0;i<fn;i++)
		if (Segments_t[i]>0)
		   for (int j=0;j<3;j++)
		   		if (Segments_v[Triangles[i][j]]==0)
		   		{
				    Segments_v[Triangles[i][j]]=Segments_t[i];
					sv[Segments_t[i]-1]++;
		        }
	
	for (int i=0;i<k1;i++)
	 	if (sv[i]>max_size)
		    k1=segment_cut(Segments_v,Vertices,vn,neighbours_vv,max_size,i+1,k1);
	  
	delete[] sv;
	sv=0;
	multiply(list_t,fn,-1);
	zeros(Segments_t,fn);				
	k2=segmentation_c(Segments_t,neighbours_tt,list_t,fn);
	sv=rzeros(sv,k2);
	//cout<<k2<<"\n";	
	for (int i=0;i<fn;i++)
		if (Segments_t[i]>0)
		   for (int j=0;j<3;j++)
		   		if (Segments_v[Triangles[i][j]]==0)
		   		{
				    Segments_v[Triangles[i][j]]=Segments_t[i]+k1;
				    sv[Segments_t[i]-1]++;
		        }
	  
	for (int i=0;i<k2;i++)
	 	if (sv[i]>max_size)
	 	{
		    k2=segment_cut(Segments_v,Vertices,vn,neighbours_vv,max_size,i+1+k1,k1+k2);
			k2-=k1;
	    }
					
	k=(int)((k1+k2)/2);
	
	
	
	for (int i=0;i<vn;i++)
		if (Segments_v[i]==0)
		{
		   	list_v[i]=1;
			vn_1++;
		}		   	
	
  Distances=rzeros(Distances,vn_1,vn_1);  
  Segments_1=rzeros(Segments_1,vn_1);
  map=rzeros(map,vn_1);
  Vertices_1=rzeros(Vertices_1,vn_1,3);
  map_1=rzeros(map_1,vn);  
  
  mapping(map,map_1,Vertices_1,vn_1,Vertices, list_v,vn);
  neighbours_vv_1=neighbours(neighbours_vv,vn,map_1,vn_1);
  
  //%перенумеровывает несегментированные вершины и записывает соответствие в
  distance_piram(Distances,Vertices_1,neighbours_vv_1,vn_1);
  //print_f(Distances,vn_1,vn_1,"Distances.txt");
  k=k_means_1(Segments_1,k,Vertices_1,Distances,vn_1);
  cout<<k1<<"\n";
  cout<<k2<<"\n";
  cout<<k<<"\n";
  for (int i=0;i<vn_1;i++)
      if (Segments_1[i]>0)	
  	  Segments_v[map[i]]=Segments_1[i]+k1+k2;
	
	delete[]Segments_t;
	delete[]Segments_1;
	delete[]list_v;
	delete[]map;
	delete[]map_1;
	delete[] list_t;
	delete[] sv;
	
	destroy(neighbours_vv_1,vn_1);
		
	destroy(Distances,vn_1);
	
	destroy(Vertices_1,vn_1);
	
 	
 	numbers[0]=k1;
 	numbers[1]=k1+k2;
 	numbers[2]=k1+k2+k;
 	return;
}
