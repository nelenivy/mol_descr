function Bmatrix=Bmarker(CMatrix)
% na vhode matrica svyaznosti
% vozvrawaet spisok b-markerov
n=size(CMatrix);
Bmatrix=[];
symb='s';

for I=1:n(1)
    for J=1:n(1)
        
        if ((CMatrix(I, J)==2) && symb=='d') % 2 dvoynie svyazi
            symb='w';
            continue;
        end
        
        if (CMatrix(I, J)==2)
            symb='d';
            continue;
        end
        
        if (CMatrix(I, J)==3) % troynaya svyaz
            symb='t';
            break;
        end
        
    end
     Bmatrix=cat(1, Bmatrix, symb);
     symb='s';
end

% s- vse sv9zi odinarnie
% d - est dvoynaya sv9z
% t- est troynaya sv9z
% a- aromaticheskaya
% w- 2 dvoynie sv9zi