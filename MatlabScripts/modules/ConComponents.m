function count=ConComponents(CMatrix)
% retun connected components quantity
str='';
str2={};
str3={};
Components={};
n=size(CMatrix);
% Components soderjit komponenty sv9znosti kak stroki

%---vydelyaetsya komponenta, soderjawaya perviy atom
% vse atomy, ne sv9zannie s pervim, zanos9tsya v str2
for I=1:n(1)
    a=TDistance(CMatrix, 1, I);
    if ((a>0) || (I==1))
       str=[str, num2str(I)]; 
    else
        str2=[str2, I];
    end
   
end
Components=[Components, str];
str='';

%---vydelyaetsya komponenta, soderjawaya perviy atom v str2
% i t.d. poka v str2 est atomy
while (size(str2)>0)
    n=size(str2);
    for I=1:n(2)
        a=TDistance(CMatrix, str2{1}, str2{I});
       if ((a>0) || (str2{I}==str2{1})) 
            str=[str, num2str(str2{I})];
        else
            str3=[str3, I];
        end
    end
    str2=str3;
    str3={};
    Components=[Components, str];
    str='';

end

n=size(Components);
count=n(2);