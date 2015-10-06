function show_surface_diff(triangles_file, prop_file_1,prop_file_2)

[vertice_matrix1, face_matrix1, prop_matrix1] = read_surface(triangles_file, prop_file_1, 1);
[vertice_matrix2, face_matrix2, prop_matrix2] = read_surface(triangles_file, prop_file_2, 1);

for i = 1:length(prop_matrix1)
    diff = prop_matrix1(i) - prop_matrix2(i);
    if (diff < 0.0)        
        prop_diff(i) = -1.0;
    elseif (diff == 0.0) 
        prop_diff(i) = 0.0;
    else  
        prop_diff(i) = 1.0;
    end
end

centers_matrix = {};
    show_surface(vertice_matrix1, face_matrix1, prop_diff, centers_matrix, 1, -1.0, 1.0);
    
end