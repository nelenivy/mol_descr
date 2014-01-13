%data �������� �� ����� result.txt ������������ � ������� ����������-������� ������������
%y �������� �� ����� y.txt ������������ � ������� ���������� - ������ ����������
% M ������� 
%1'euclidean' � Euclidean distance 
%2'cityblock' � Sum of absolute differences 
%3'cosine' � One minus the cosine of the included angle between points (treated as vectors)
%4'correlation' � One minus the sample correlation between points (treated as sequences of values)
%5'hamming' � Percentage of bits that differ (only suitable for binary data)
%MAX - ���-�� ������������
%N - ����� ������������� �����������
%kol - ���-�� ������������, ������� ����� ��������
function [res] = estumate_v4(M,MAX,n_molecules,kol,koef,l_granutcua,r_granutcua,glybuna,kol_coced,sampl,vector,y,Rad)
num=0;
t=1;
max=-1000000000;
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
        R_2(i)= 0;
        otkaz(i)=0;
    end
    for i=1:kol_coced
        y_pr=0;
        y_sr=0;
        for k=1:n_molecules
            if(Class(k,i)==1000000)
                otkaz(i)=otkaz(i)+1;
                continue;
            end
            y_sr =y_sr+y(k);
        end
        y_sr=y_sr/(n_molecules-otkaz(i));
        for k=1:n_molecules
            if(Class(k,i)==1000000)
                continue;
            end
            R_2(i)= R_2(i)+ (Class(k,i)-y(k))*(Class(k,i)-y(k));
            y_pr =y_pr+(y_sr-y(k))*(y_sr-y(k));
        end
        R_2(i)= 1-(R_2(i)/y_pr);
    end
    R_2
    for i1=1:kol_coced
        neopr=otkaz(i1);
        Da=0;
       % for i=1:n_molecules
        %    if Class(i,i1)>l_granutcua&&Class(i,i1)<r_granutcua&&(Class(i,i1)+y(i)>0||Class(i,i1)+y(i)<0)
         %       Da=Da+1;
          %  end
           %if Class(i,i1)+y(i)==0
            %    neopr=neopr+1;
         %   end
       % end
     %   if koef==2
      %      if(n_molecules-neopr==0)
       %         count=0;
        %    else
         %   count=Da/(n_molecules-neopr)*100;
          %  end
   %     else
    %    count=Da+koef*neopr;
     %   end
        count=R_2(i1);
        if count>=max
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
    max=-1000000000000;
end
res=result;
