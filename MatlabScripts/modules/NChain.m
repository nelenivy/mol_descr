function chains=NChain(file,k)

if isempty(file) || isempty(k)
    return
end
chains=[];
chainsNew=[];

CMatrix=ConnectMatrix(file);

n=size(CMatrix);
% 2 atoma v cepochke
for I=1:n(1)
    for J=I+1:n(2)
        if (CMatrix(I,J)~=0)
            temp(1)=I;
            temp(2)=J;
            chains=[chains; temp];
            temp(1)=J;
            temp(2)=I;
            chains=[chains; temp];
        end
    end
end


for length=2:k-1
    m=size(chains);
    for I=1:m(1)
        for J=1:n(2)
             if ((CMatrix(chains(I, length),J)~=0) && all(J~=chains(I, 1:length))) % iskluchaem stroki vida a-b-a
                 chains(I, length+1)=J;
                 chainsNew=[chainsNew; chains(I,:)];
             end
        end
    end
    chains=chainsNew;
    chainsNew=[];
end

chains_numb=ReversEqual(chains);

% сохранить результат
ind=find(file=='\',1,'last');
warning off
mkdir(strcat(file(1:ind),'linear_fragments\level1\',num2str(k)));
warning on
save(strcat(file(1:ind),'linear_fragments\level1\',num2str(k),'\chains_numb.mat'),'chains_numb');