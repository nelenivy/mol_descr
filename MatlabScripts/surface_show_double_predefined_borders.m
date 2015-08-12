function  surface_show_double_predefined_borders(triangles_file, prop_file, centers_file, double_prop, use_borders, min_val, max_val)

fid=fopen(triangles_file,'r');

input_buf = fscanf(fid,'%*c %f %*s %f %*s %f %*s',[1, 3]);  
fgets(fid);                   
face_ind = 1;	
while(~(isempty(input_buf)))    
        face_matrix(face_ind,1) = input_buf(1) + 1;
        face_matrix(face_ind,2) = input_buf(2) + 1;
        face_matrix(face_ind,3) = input_buf(3) + 1; 
        input_buf = fscanf(fid,'%*c %f %*s %f %*s %f %*s',[1, 3]); 
        fgets(fid); 
        face_ind = face_ind + 1;
end
fclose(fid);

vert_ind = 1;
fid=fopen(prop_file,'r');
prop_format = '%*c %f %*s %f %*s %f %*s %d %*s';
if (double_prop)
   prop_format = '%*c %f %*s %f %*s %f %*s %f %*s';
end
input_buf = fscanf(fid,prop_format,[1, 4]); 
fgets(fid);   
max_segment = 0;
min_segment = 100000;
while(~(isempty(input_buf)))
        vertice_matrix(vert_ind,1) = input_buf(1);
        vertice_matrix(vert_ind,2) = input_buf(2);
        vertice_matrix(vert_ind,3) = input_buf(3);
        prop_matrix(vert_ind) = input_buf(4);
        
        if (prop_matrix(vert_ind) > max_segment)
           max_segment =  prop_matrix(vert_ind);
        end
        
        if (prop_matrix(vert_ind) < min_segment)
           min_segment =  prop_matrix(vert_ind);
        end
        
        input_buf=fscanf(fid,prop_format,[1, 4]); 
        fgets(fid);
        vert_ind = vert_ind + 1;
end

fclose(fid);
max_segment
min_segment

if (use_borders ~= 1)
    max_val = max_segment;
    min_val = min_segment;
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

fid=fopen(prop_file,'r');
input_buf=fscanf(fid,'%*c %f %*s %f %*s %f %*s %f %*s',[1, 4]); 
fgets(fid);   
type_ind = 1;
centers_num = length(centers_matrix(:, 1))

while(~(isempty(input_buf))) 
        colour_matrix(type_ind,1) = 0;
        colour_matrix(type_ind,2) = 0;
        colour_matrix(type_ind,3) = 0;
        
        found = 0;
        
        for i = 1 : centers_num
            if (vertice_matrix(type_ind, :) == centers_matrix(i, :))
                %%vertice_matrix(type_ind, :)
                %%centers_matrix(i, :)
                found = 1;
                break;
            end
        end
        
        if (found == 1)
            colour_matrix(type_ind,1) = 1;
            colour_matrix(type_ind,2) = 1;
            colour_matrix(type_ind,3) = 1;
        else  
            prop = input_buf(4);
            if (double_prop == 0)
               prop = prop / max_segment;
            else
                prop = (prop-min_val) / (max_val - min_val);
                prop = min (prop, 1.0);
                prop = max (prop, 0.0);
            end  
            prop = prop * (255.0 - 85.0) + 85.0;
            [R, G, B] = HSVToRGB(prop, 240, 240);
            %curr_color = max(ceil(prop * 3.0), 1);
            %color_val = (prop - (curr_color - 1) / 3.0) * 3.0;      
            %color_val = max(color_val, 0.1);
            %colour_matrix(type_ind,curr_color) = color_val;
            colour_matrix(type_ind,1)=R;
            colour_matrix(type_ind,2)=G;
            colour_matrix(type_ind,3)=B;
        end
            
        input_buf=fscanf(fid,'%*c %f %*s %f %*s %f %*s %f %*s',[1, 4]); 
        fgets(fid);
        type_ind = type_ind + 1;
end

fclose(fid);
size(vertice_matrix)
size(face_matrix)
size(colour_matrix)
patch('Vertices',vertice_matrix,'Faces',face_matrix,'FaceVertexCData',colour_matrix,'FaceColor','interp','FaceAlpha',1);
 view(3);
