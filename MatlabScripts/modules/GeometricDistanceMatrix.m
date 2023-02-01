function GDistMatrix=GeometricDistanceMatrix(file)
% vozvrawaet matricu geometricheskih rasstoyaniy

if (file=='0')
   return
end
if (isempty(file))
    return
end

[Coordinates, Bound, AtomName]=molfile2matrixes(file);
n=size(Coordinates);
GDistMatrix=zeros(n(1));
for I=1:n(1)
    for J=I+1:n(1)
        distx=(Coordinates(I,1)-Coordinates(J,1))^2;
        disty=(Coordinates(I,2)-Coordinates(J,2))^2;
        distz=(Coordinates(I,3)-Coordinates(J,3))^2;
        GDistMatrix(I,J)=sqrt(distx+disty+distz);
        GDistMatrix(J,I)=GDistMatrix(I,J);
    end
end
