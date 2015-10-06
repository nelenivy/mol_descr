function show_surface(vertice_matrix, face_matrix, prop_matrix, centers_matrix, double_prop, min_val, max_val)

if (~isempty(centers_matrix))
    centers_num = length(centers_matrix(:, 1));
else
    centers_num = 0;
end

vert_number = length(prop_matrix);
for type_ind = 1 : vert_number
        colour_matrix(type_ind,1) = 0;
        colour_matrix(type_ind,2) = 0;
        colour_matrix(type_ind,3) = 0;
        
        found = 0;
        
        for i = 1 : centers_num
            if (vertice_matrix(type_ind, :) == centers_matrix(i, :))
                found = 1;
                break;
            end
        end
        prop = (prop_matrix(type_ind));
        if (found == 1)
            colour_matrix(type_ind,1) = 1;
            colour_matrix(type_ind,2) = 1;
            colour_matrix(type_ind,3) = 1;
        elseif (0)%prop < 0.5 && prop> -0.5)
               colour_matrix(type_ind,1)=1;
            colour_matrix(type_ind,2)=1;
            colour_matrix(type_ind,3)=0; 
        else             
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
end

size(vertice_matrix)
size(face_matrix)
size(colour_matrix)
patch('Vertices',vertice_matrix,'Faces',face_matrix,'FaceVertexCData',colour_matrix,'FaceColor','interp','FaceAlpha',1);
 view(3);

end