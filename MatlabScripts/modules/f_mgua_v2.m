% MGUA process
% for matrix of properties 'xmgua' and row vector 'ymgua'

% Required values in the workspace:
% q_single - mgua depth for the first selection
% q_mult - mgua depth for the other selections
% vartol - variance threshold
% corrtol_single - correlation threshold for the first selection
% corrtol_mult - correlation threshold for the other selections
% no_best_models - the number of the best models saved

function [result, used, coefficients, R2, informative] = f_mgua_v2(xmgua, ymgua, q_single, q_mult, vartol, corrtol_single, corrtol_mult, no_best_models)

% OUTCOME:
% result - vectors 'b' for 'no_best_models' best models (no_best_models x (m + 1))
% use = (result == 0) (no_best_models x (m + 1))
% coefficients - coefficients of one best model (without zeros)
% R2 - for one best model
% informative - numbers of features involved in the best models


s = size (xmgua);
m = s(2); % - number of objects' properties
n = s(1); % - number of objects

if (length(unique(ymgua)) == 2)
    isClassification = 1;
else
    isClassification = 0;
end


% NORMALISATION
stdx = std (xmgua, 1); % - standard deviation
meanx = mean (xmgua); % - mean value
stdx(stdx == 0) = 1;

xmgua_original = xmgua;
ymgua_original = ymgua;

for i = m : (-1) : 1
    xmgua(:, i) = (xmgua (:, i) - meanx (i))/ stdx (i);
end

stdy = std (ymgua, 1);
meany = mean (ymgua);
ymgua = (ymgua - meany) / stdy;

i=0;
j=0;
k=0;

single = [];
multiple = [];


% 'single' construction
display('single')

while i < m
    i = i + 1;
    x1 = [ones(n,1) xmgua(1:n, i)];
    k = k + 1;
    [a,bint1,r1,rint1,stats1] = regress(ymgua', x1);
    if(isequal(a,0) && isequal(bint1,0) && isequal(r1,0) && isequal(rint1,0) && isequal(stats1,0))
        result=[0]; used=[]; coefficients=[]; R2=[]; informative=[];
        return;
    end
    y1 = x1 * a;

    b = sum(y1' .* ymgua < 0);
%   OLD VERSION
%     b = length(y1);
% 
%     for l = 1 : length(y1)
%         b = b - (y1(l) * ymgua(l) > 0);
%     end

    b2 = ((y1 - ymgua')' * (y1 - ymgua'));
    new_col = [a; b; i; y1; b2];
    single = [single new_col];
end   

if isClassification == 1
    single = (sortrows (single', [3 (n + 5)]))';
else
    single = (sortrows (single', [n + 5]))';
end

xfilter = single (5 : (n + 4), :);
[xfilter, filtered] = f_filter (xfilter, vartol, corrtol_single);

single (:, filtered) = [];
% for i = 1 : filtnum
%     single (:, filtered (i)) = [];
% end
if (size(single, 2) > q_single) 
    single = single (:, 1 : q_single);
end


% first 'multiple' construction

display('first multiple')

i = 0;
j = 0;
k = 0;

s = 1;
loc_multiple = [];
while i < size(single, 2)
    i = i + 1;
    j = i;
    while j <= size (single, 2)
        x1 = [ones(n,1) single(5 : (n + 4), [i j])];
        k = k + 1;
        a = regress(ymgua', x1);
        y1 = x1 * a;
        
        b = sum(y1' .* ymgua < 0);
        
%         b = length(y1);
%         for l = 1 : length(y1)
%             b = b - (y1(l) * ymgua(l) > 0);
%         end
        b2 = ((y1 - ymgua')' * (y1 - ymgua'));
        new_col = [a; b; i; j; y1; b2];
        loc_multiple(1 : (n + 7), k) = new_col;
        j = j + 1;
    end
end

if isClassification == 1
    loc_multiple (1 : (n + 7), 1 : k)= (sortrows ((loc_multiple(1 : (n + 7),1 : k))' , [4 (n + 7)]))';
else
    loc_multiple (1 : (n + 7), 1 : k)= (sortrows ((loc_multiple(1 : (n + 7),1 : k))' , [n + 7]))';
end
xfilter = loc_multiple (7 : (n + 6), :);

[xfilter, filtered] = f_filter (xfilter, vartol, corrtol_mult);

loc_multiple(:, filtered) = [];
% for i = 1 : filtnum
%     loc_multiple(:, filtered (i)) = [];
% end
multiple(1 : (n + 7), 1 : min(q_mult, size(loc_multiple, 2)), s)= loc_multiple(1: (n + 7), 1 : min(q_mult, size(loc_multiple, 2)));
multiple_size = min(q_mult, size(loc_multiple, 2));
s = s + 1;
   

% following 'multiple' construction

while s < (n / 5) % - number of properties used in mgua
%while s < min (n/5, 3)
    display(s)
    
    i=0; % 
    j=0; % 
    k=0; % - current equation number
    loc_multiple = [];
    
    for i = 1 : size(single, 2)
        for j = 1 : multiple_size
            x1 = [ones(n,1) single(5 : (n + 4), i) multiple(7 : (n + 6), j,  s - 1)];
            k = k + 1;
            a = regress(ymgua', x1);
            y1 = x1 * a;
            b = sum(y1' .* ymgua < 0);
            
%             b = length(y1);
%             for l = 1 : length(y1)
%                 b = b - (y1(l) * ymgua(l) > 0);
%             end
            b2 = ((y1 - ymgua')' * (y1 - ymgua'));
                       
            new_col = [a; b; i; j; y1; b2];
            loc_multiple(1 : (n + 7), k) = new_col;
        end
    end
        
    if isClassification == 1
        loc_multiple(1 : (n + 7), 1 : k)= (sortrows ((loc_multiple(1 : (n + 7),1 : k))', [4 (n + 7)]))';
    else
        loc_multiple(1 : (n + 7), 1 : k)= (sortrows ((loc_multiple(1 : (n + 7),1 : k))', [n + 7]))';
    end
    
    xfilter = loc_multiple(7 : (n + 6), :);
    
    [xfilter, filtered] = f_filter (xfilter, vartol, corrtol_mult);
    
    loc_multiple(:, filtered) = [];
%     for i = 1 : filtnum
%         loc_multiple(:, filtered (i)) = [];
%     end
    multiple(1 : (n + 7), 1 : min(q_mult, size(loc_multiple, 2)), s)= loc_multiple(1: (n + 7), 1 : min(q_mult, size(loc_multiple, 2)));
    multiple_size = min(q_mult, size(loc_multiple, 2));

    s = s + 1;
end
        


% return trace

no_best_models = min(no_best_models, multiple_size);

result = zeros (m + 1, no_best_models);
used = zeros (m, no_best_models);

for start_index = 1 : no_best_models
    
    %res = zeros (1, size(single, 2) + 1);
    res = zeros(1, m + 1);
    index = start_index;
    
    c = 1;
    for i = (s - 1) : -1 : 1
        res (size(single, 2)+1) = res (size(single, 2)+1) + c * multiple(1, index, i); % свободный член
        res (multiple(5, index, i)) = res(multiple(5, index, i)) + c * multiple(2, index, i);
        c = c * multiple(3, index, i);
        if (i~=1)
            index = multiple(6, index, i);
        end
    end

    res(multiple(6, index, 1))= res(multiple(6, index, 1)) + c; % multiple! not single coefficients


    result (m + 1, start_index) = res (size(single, 2) + 1);
    for i = 1 : size(single, 2)
        result(m + 1, start_index) = result(m + 1, start_index) + res(i) * single(1, i);
        result(single (4, i), start_index) = result (single (4, i), start_index) + res (i) * single(2, i);
    end

    
    used(result(1:m, :) ~= 0) = 1;

    % for i = 1 : m
    %     if (result(i) ~= 0)
    %         used(i) = 1;
    %     else
    %         used(i) = 0;
    %     end
    % end

    %result
    %used

    % checking
    error = ([xmgua ones(n, 1)] * result(:, start_index) - ymgua')' * ([xmgua ones(n, 1)] * result(:, start_index) - ymgua') - multiple (n + 7, start_index, size(multiple, 3))
end


coefficients = [];

i = 0;
for j = 1 : size(xmgua, 2)
    if (used(j, 1) == 1)
        i = i + 1;
        coefficients(i)= result(j, 1); %coefficient
    end

end

coefficients(i + 1) = result(size(xmgua, 2) + 1, 1);

% correct coefficients because of scaling
coefficients(i + 1) = coefficients(i + 1) * stdy + meany;
for j  = 1 : i
    coefficients(j) = coefficients(j)* stdy / stdx (j);
    coefficients(i + 1) = coefficients(i + 1) - (coefficients(j) * meanx(j));
end

xmgua = xmgua_original;
ymgua = ymgua_original;

xmgua(:, used(:, 1) == 0) = [];
xmgua = [xmgua ones(n, 1)];

control = 0;
if isClassification == 1
    for i = 1 : size (ymgua, 2)
        x_reduced = xmgua;
        y_reduced = ymgua;
        x_reduced (i, :) = [];
        y_reduced (i) = [];
%         if isempty(find(used ~=0, 1))
%             control = control + (coefficients(clustnum, 1) * ymgua(i) >= 0);
%         else
        control  = control + ((xmgua(i, :) * regress (y_reduced', x_reduced)) * ymgua(i) >= 0);
%        end
    end
    R2 = control / size (ymgua, 2);

else

    for i = 1 : size (ymgua, 2)
        x_reduced = xmgua;
        y_reduced = ymgua;
        x_reduced (i, :) = [];
        y_reduced (i) = [];
        control = control + (xmgua(i, :) * regress (y_reduced', x_reduced) - ymgua(i)) ^ 2;
    end


    R2 = 1 - (control / ((ymgua - mean(ymgua)) * (ymgua - mean(ymgua))'));
end

informative = find(max(used, [], 2));
