function show_scale_extremums(triangles_file, prop_file,prop_file_prev, prop_file_next)

[vertice_matrix, face_matrix, prop_matrix] = read_surface(triangles_file, prop_file, 1);
[vertice_matrix1, face_matrix1, prop_matrix_prev] = read_surface(triangles_file, prop_file_prev, 1);
[vertice_matrix2, face_matrix2, prop_matrix_next] = read_surface(triangles_file, prop_file_next, 1);

for i = 1:length(prop_matrix_prev)
    diff_prev = prop_matrix(i) - prop_matrix_prev(i);
    diff_next = prop_matrix(i) - prop_matrix_next(i);
    if (diff_prev * diff_next >= 0.0)        
        scale_extr(i) = 1.0;
    else
        scale_extr(i) = -1.0;
    end
end

centers_matrix = {};
    show_surface(vertice_matrix1, face_matrix1, scale_extr, centers_matrix, 1, -1.0, 1.0);
    
end