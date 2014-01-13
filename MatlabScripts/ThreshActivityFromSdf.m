function ThreshActivityFromSdf(sdf, output_file, activity_string, threshold)
activity_string;
fin=fopen(sdf,'r');
fout=fopen(output_file,'w');
line=fgets(fin);
size(line);
found = 0;

while(line~=-1)  
    if(found == 1)        
        str2double(line)
		num = str2double(line);
		
		if (num < threshold)
			fprintf(fout,'-1 \n');
		else
			fprintf(fout,'1 \n');
		end
		
        found = 0;
    end
	
    if(findstr(activity_string,line) == 1)
        found = 1;
    end
    
    line=fgets(fin);
end

fclose(fin);
fclose(fout);