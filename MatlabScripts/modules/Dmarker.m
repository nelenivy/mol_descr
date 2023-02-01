function Dmatrix=Dmarker(CMatrix)
% d-marker - kolichestvo sv9zey atoma
n=size(CMatrix);
Dmatrix=[];
count=0;

for I=1:n(1)
    for J=1:n(1)
        if (CMatrix(I, J)~=0)
            count=count+1;
        end
    end
     Dmatrix=[Dmatrix; num2str(count)];
     count=0;
end