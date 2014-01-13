function DiskrActivityFromSdf(sdf, output_file, activity_string, no_activity_string)
fin=fopen(sdf,'r');
fout=fopen(output_file,'w');
line=fgets(fin);
size(line);
while(line~=-1)
 %l=size(line)   
    %if(i==1)
        %fprintf(fout,line);
        %i=0;
   % end
    
    %if(l(2)==14)
        if(findstr(activity_string, line) == 1)
            fprintf(fout,'1 \n');
        elseif(findstr(no_activity_string, line) == 1)
            fprintf(fout,'-1 \n');
        end
    %end
    line=fgets(fin);
end
fclose(fin);
fclose(fout);