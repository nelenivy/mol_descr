function [vertice_matrix, face_matrix, prop_matrix] = read_surface(triangles_file, prop_file, double_prop, sphere_prop)

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

while(~(isempty(input_buf)))
        vertice_matrix(vert_ind,1) = input_buf(1);
        vertice_matrix(vert_ind,2) = input_buf(2);
        vertice_matrix(vert_ind,3) = input_buf(3);
        prop_matrix(vert_ind) = input_buf(4);      
        
        input_buf=fscanf(fid,prop_format,[1, 4]); 
        fgets(fid);
        vert_ind = vert_ind + 1;
end

fclose(fid);

end
