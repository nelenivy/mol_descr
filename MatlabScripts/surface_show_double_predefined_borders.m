function  surface_show_double_predefined_borders(triangles_file, prop_file, centers_file, double_prop, use_borders, min_val, max_val, use_dev)

[vertice_matrix, face_matrix, prop_matrix] = read_surface(triangles_file, prop_file, double_prop);

%calculate mean and deviation of prop_matrix
max_segment = 0;
min_segment = 100000;
mean = 0.0;
for i = 1:length(prop_matrix)
   mean  = mean + prop_matrix(i);
   if (prop_matrix(i) > max_segment)
       max_segment =  prop_matrix(i);
    end

    if (prop_matrix(i) < min_segment)
       min_segment =  prop_matrix(i);
    end
end
mean = mean / length(prop_matrix);
dev = 0.0;
for i = 1:length(prop_matrix)
   dev = dev + (prop_matrix(i) - mean)*(prop_matrix(i) - mean);
end
dev = dev / length(prop_matrix)
mean = -2.0;

if (~use_borders) 
    min_val = min_segment
    max_val = max_segment
end

if (use_dev)
    min_val = mean - 4.0 * sqrt(dev)
    max_val = mean + 4.0 * sqrt(dev)
end

fid=fopen(centers_file,'r');

input_buf = fscanf(fid,'%*c %f %*s %f %*s %f %*s %d %*s',[1, 4]); 
fgets(fid);
center_ind = 1;

while(~(isempty(input_buf)))
        centers_matrix(center_ind,1) = input_buf(1);
        centers_matrix(center_ind,2) = input_buf(2);
        centers_matrix(center_ind,3) = input_buf(3);        
           
        input_buf=fscanf(fid,'%*c %f %*s %f %*s %f %*s %d %*s',[1, 4]); 
        fgets(fid);
        center_ind = center_ind + 1;
end
fclose(fid);

show_surface(vertice_matrix, face_matrix, prop_matrix, centers_matrix, double_prop, min_val, max_val);
