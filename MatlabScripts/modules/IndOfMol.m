function ind=IndOfMol(Molecula, Matrix)
ind=0;
for I=1:size(Matrix,1)
    if (Molecula==Matrix(I,:))
        ind=I;
        break;
    end
end