function res=binary_search(Array, Str)

% binarnym poiskom iwet stroku Str v massive Array
% esli stroka est v massive, vozvrawet nomer stroki
% net -  0
res=0;
ArraySize=size(Array);
if (ArraySize==0)
    return;
end
left=1;
right=ArraySize(1);
ind=ceil((right-left)/2);
% sravnit s pervim i poslednim
if (strcmp(Array{left}, Str))
    res=1;
    return;
end
if (strcmp(Array{right}, Str))
    res=right;
    return;
end
if (ind==0)
    return;
end

while (ind ~=1 && ind~=ArraySize(1))
    if (strcmp(Array{ind}, Str))
        res=ind;
        return;
    else
        
        StrSize=size(Str);
        ArrayStr=size(Array{ind});
        StrSizeMin=min(StrSize,ArrayStr);
        for I=1:StrSizeMin(2)
            if (Str(I)==Array{ind}(I))
                res=ind;
            end
            if (Str(I)<Array{ind}(I))
                right=ind; % pervaya chast massiva
                ind=left+ceil((right-left)/2);
                if (right==left+1)
                    ind=ceil((right-left)/2);
                end
                res=0;
                break;
            end
            if (Str(I)>Array{ind}(I))
                left=ind; % vtoraya chast massiva
                ind=ind+ceil((right-left)/2);
                res=0;
                break;
            end
        end
        % esli Str==Array{ind}
        if (res==ind && (StrSize(2)==ArrayStr(2)))
            return;
        end
        if (res==ind && StrSize(2)>ArrayStr(2))
                left=ind; % vtoraya chast massiva
                ind=ind+ceil((right-left)/2);
                res=0;                
        end
        if (res==ind &&StrSize(2)<ArrayStr(2))
                left=ind; % vtoraya chast massiva
                ind=ind+ceil((right-left)/2);
                res=0;
        end
        
    end
end