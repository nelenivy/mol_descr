%data вводитс€ из фаила result.txt наход€щегос€ в текущей дериктории-матрица дескрипторов
%y вводитс€ из фаила y.txt наход€щегос€ в текущей дериктории - вектор активности
% M метрика 
%1'euclidean' Ч Euclidean distance 
%2'cityblock' Ч Sum of absolute differences 
%3'cosine' Ч One minus the cosine of the included angle between points (treated as vectors)
%4'correlation' Ч One minus the sample correlation between points (treated as sequences of values)
%5'hamming' Ч Percentage of bits that differ (only suitable for binary data)
%MAX - кол-во дескрипторов
%N - номер фиксированого дескриптора
%kol - кол-во дескрипторов, которое будет отобрано
function [res] = estumate_v4(M,MAX,n_molecules,kol,koef,l_granutcua,r_granutcua,glybuna,kol_coced,sampl,vector,y,Rad)
num=0;
t=1;
max=0;
J=0;
k1=0;
A=size(vector)-4;
N=A(1);
%M='euclidean';
for i=1:n_molecules
    for j=1:kol_coced
        Class(i,j)=0;
    end
end
for i=1:N+5
    for j=1:kol
        result(i,j)=0;
    end
end
DA=0;
for j=1:MAX
    Da=0;
    Neopr=0;
    sampl1=sampl(j,:);
    for j1=1:N
        sampl1=[sampl1;sampl(vector(j1+4),:)]; 
    end
    if Rad~=0
        raduys=opr_raduys(sampl1,Rad);
    else
        raduys=10000000000000;
    end
    training =sampl1(:,2:n_molecules)';
    group = [y(2:n_molecules)];
    Molekula=sampl1(:,1)';
    Class(1,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    for i=2:n_molecules-1
        training = [sampl1(:,1:i-1),sampl1(:,i+1:n_molecules)]';
        group = [y(1:i-1);y(i+1:n_molecules)];
        Molekula=sampl1(:,i)';
        Class(i,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    end    
    training = sampl1(:,1:n_molecules-1)';
    group = y(1:n_molecules-1);
    Molekula=sampl1(:,n_molecules)';
    Class(n_molecules,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    for i=1:kol_coced
        Class(:,i) =Class(:,i)-y;
    end
    for i1=1:kol_coced
        neopr=0;
        Da=0;
        for i=1:n_molecules
            if Class(i,i1)>l_granutcua&&Class(i,i1)<r_granutcua&&(Class(i,i1)+y(i)>0||Class(i,i1)+y(i)<0)
                Da=Da+1;
            end
            if Class(i,i1)+y(i)==0
                neopr=neopr+1;
            end
        end
        if koef==2
            if(n_molecules-neopr==0)
                count=0;
            else
            count=Da/(n_molecules-neopr)*100;
            end
        else
        count=Da+koef*neopr;
        end
        if count>max
            max=count;
            DA=Da;
            Neopr=neopr;
            k1=i1;
            J=j;
        end
       % Da
    end
    if num < kol
        num=num+1;
        result(1,num)=max;
        result(2,num)=DA;
        result(3,num)=Neopr;
        result(4,num)=k1;
        for j1=1:N
            result(4+j1,num)=vector(4+j1);
        end
        result(5+N,num)=J;        
    else
        I=1;
        lok_min=result(1,1);
        for i=2:kol
            if lok_min>result(1,i)
                I=i;
                lok_min=result(1,i);
            end
        end
        if result(1,I) < max
            result(1,I)=max;
            result(2,I)=DA;
            result(3,I)=Neopr;
            result(4,I)=k1;
            result(5+N,I)=J;
        end         
    end
    max=0;
end
res=result;
