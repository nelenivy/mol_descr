function ConnectMatrix=ConnectMatrix(file)
% return connectivity matrix
if (file=='0')
   return
end
if isempty(file)
    return
end

[Coordinates, Bound, Atoms]=molfile2matrixes(file);
n=size(Bound);
m=size(Atoms);
ConnectMatrix=zeros(m(1));
for I=1:n(1)
    ConnectMatrix(Bound(I,1),Bound(I,2))=Bound(I,3);
    ConnectMatrix(Bound(I,2), Bound(I,1))=Bound(I,3);
end
