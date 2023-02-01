% filter to skip constant and almost constant columns and columns with big
% value of correlation
% CHANGE There is no filtering of columns with big value of correlation

% Required values in the workspace:
% vartol - variance threshold
% corrtol - correlation threshold


function [xfilter, filtered] = f_filter (xfilter, vartol, corrtol)

loc_size = size (xfilter);
i = loc_size (2);
filtnum = 0;
filtered = [];
% if issingle == 1
%     corrtol = corrtol_single;
% else
%     corrtol = corrtol_mult;
% end

while i > 0
    if (var (xfilter (:, i)) < vartol)
        xfilter (:, i) = [];
        filtnum = filtnum + 1;
        filtered (filtnum) = i;
    else
        for j = 1 : i - 1
            if (corr (xfilter (:, j), xfilter(:, i)) > corrtol)
                xfilter (:, i) = [];
                filtnum = filtnum + 1;
                filtered (filtnum) = i;
                break;
            end
        end
    end
    i = i - 1;
end