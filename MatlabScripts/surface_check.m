function  surface_check(triangles_file, segments_file, types_file, centers_file)

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
fid=fopen(segments_file,'r');
input_buf = fscanf(fid,'%*c %f %*s %f %*s %f %*s %d %*s',[1, 4]); 
fgets(fid);   
max_segment = 0;

while(~(isempty(input_buf)))
        vertice_matrix(vert_ind,1) = input_buf(1);
        vertice_matrix(vert_ind,2) = input_buf(2);
        vertice_matrix(vert_ind,3) = input_buf(3);
        segments_matrix(vert_ind) = input_buf(4);
        
        if (segments_matrix(vert_ind) > max_segment)
           max_segment =  segments_matrix(vert_ind);
        end
        
        input_buf=fscanf(fid,'%*c %f %*s %f %*s %f %*s %d %*s',[1, 4]); 
        fgets(fid);
        vert_ind = vert_ind + 1;
end

fclose(fid);
max_segment

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

fid=fopen(types_file,'r');
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
            curr_type = fix((input_buf(4) - 1) / 3) + 1
            curr_type1 = mod(input_buf(4) - 1, 3) + 1;
            
            if (curr_type > 3)
                colour_matrix(type_ind,curr_type) = segments_matrix(type_ind) / max_segment;
            end
            colour_matrix(type_ind,curr_type1) = segments_matrix(type_ind) / max_segment;
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
