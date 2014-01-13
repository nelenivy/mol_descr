function ActivityFromSdf(sdf)
fin=fopen(sdf,'r');
fout=fopen('y.txt','w');
i=0;
line=fgets(fin);
size(line);
while(line~=-1)
 l=size(line);   
    if(i==1)
        fprintf(fout,line);
        i=0;
    end
    
   
    if(findstr('>  <Tg (K)exp>',line)==1)
        i=1
    end
    
    line=fgets(fin);
end
fclose(fin);
fclose(fout);
        
        