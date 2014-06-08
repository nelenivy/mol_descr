function ActivityFromSdf(sdf, output_file, activity_string)
activity_string;
fin=fopen(sdf,'r');
fout=fopen(output_file,'w');
line=fgets(fin);
size(line);
found = 0;

while(line~=-1)  
    if(found == 1)        
        str2double(line);
		num = str2double(line)		
		fprintf(fout,'%f \r\n', num)	
        found = 0;
    end
	
    if(findstr(activity_string,line) == 1)
        found = 1;
    end
    
    line=fgets(fin);
end

fclose(fin);
fclose(fout);