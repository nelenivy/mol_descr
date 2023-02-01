function ChainCode(FileName, k, format)
%vozvrawaet spisok kodirovannih cepo4ek dliny k v zadannom formate

CodeMatr=[];temp=[];Dmatr=[];Bmatr=[];Rmatr=[];NewFileName=[];
[Coord, Bound, AtomName]=molfile2matrixes(FileName);
CMatrix=ConnectMatrix(FileName);
m=size(CMatrix);

ind=find(FileName=='\',1,'last');
if (fopen(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\chains_numb.mat'),'r')==-1)
    NChain(FileName, k);
    fclose('all');
end
load(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\chains_numb.mat'));
n=size(chains_numb);

%--------------- razbor formata -----------

    if (format(1)~='_')
        Dmatr=Dmarker(CMatrix);
    else
        Dmatr(1:m(1))='_';
    end
    if (format(2)~='_')
        Bmatr=Bmarker(CMatrix);
    else
        Bmatr(1:m(1))='_';
    end
    if (format(3)~='_')
        Rmatr=Rmarker(CMatrix);
    else
        Rmatr(1:m(1))='_';
    end

%------------------------------------------
%---dlya kajdogo atoma formiruyuts9 markery
%---atom: imya, D-marker, B-marker, R-marker

for I=1:n(1)
    for J=1:k
        % a - kod 1 atoma
        a=AtomName(chains_numb(I,J),:); %imya atoma
        a=cat(2, a, Dmatr(chains_numb(I,J)));
        a=cat(2, a, Bmatr(chains_numb(I,J)));
        a=cat(2, a, Rmatr(chains_numb(I,J)));
        temp=cat(2, temp, a); % v temp sobiraetsya cepochka
    end
    CodeMatr=cat(1, CodeMatr, temp);
    temp=[];
end
chains_coded=CodeMatr;
%-------------to file -------------------

ind=find(FileName=='\',1,'last');
warning off
mkdir(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\',format));
warning on
save(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\',format,'\chains_coded.mat'),'chains_coded');