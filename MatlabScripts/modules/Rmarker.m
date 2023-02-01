function Rmatrix=Rmarker(CMatrix)

Rmatrix=[];
RmatrixNew=[];
n=size(CMatrix);
R='c';
components=ConComponents(CMatrix);
tempMatr=CMatrix;
r=0;
c=0;

% zapolnenie matricy
BondMatrix=zeros(n(1),n(1));
for I=1:n(1)
    for J=I+1:n(1)
        if ((CMatrix(I,J)~=0) && (BondMatrix(I,J)==0))
            tempMatr(I, J)=0;
            tempMatr(J, I)=0;
            newComponents=ConComponents(tempMatr);
            if (newComponents==components)
                BondMatrix(I,J)='r';
                BondMatrix(J,I)='r';
            end
            if (newComponents==components+1)
                BondMatrix(I,J)='c';
                BondMatrix(J,I)='c';
            end
        end
        tempMatr=CMatrix;
    end
end

% podschet kolichestva r i c
R='c';
for I=1:n(1)
    r=0; c=0;
    for J=1:n(1)
        if (BondMatrix(I,J)=='c')
            c=c+1;
        end
        if (BondMatrix(I,J)=='r')
            r=r+1;
        end
    end
    if ((r>=2) && (c>0))
        R='s';
    end
    if (c==0)
        R='r';
    end
    Rmatrix=cat(1, Rmatrix, R);
    R='c';
    
end
