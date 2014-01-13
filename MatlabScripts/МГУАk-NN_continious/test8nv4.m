%data вводитс€ из фаила result.txt наход€щегос€ в текущей дериктории-матрица дескрипторов
%y вводитс€ из фаила y.txt наход€щегос€ в текущей дериктории - вектор активности
% M метрика 
%1'euclidean' Ч Euclidean distance 
%2'cityblock' Ч Sum of absolute differences 
%3'cosine' Ч One minus the cosine of the included angle between points (treated as vectors)
%4'correlation' Ч One minus the sample correlation between points (treated as sequences of values)
%5'hamming' Ч Percentage of bits that differ (only suitable for binary data)
%MAX - кол-во дескрипторов
%kol - количество дескрипторов, которое будет отобрано при каждом проходе 
%результаты первого прохода записываютс€ в фаил out1.txt 
%результаты всех вторых проходов записываютс€ в фаил out2.txt
%результаты всех третьих проходов записываютс€ в фаил out3.txt и т. д.
function [] = test8nv4 (M,MAX,n_molek,kol,koef,l_granutcua,r_granutcua,glybuna,kol_coced,Rad)
%n_molecules = MAX;
n_molecules=n_molek;
%M='euclidean';
myfid =  fopen ('result.txt', 'r');%fullfile (matlabroot, 'work', 'NewFolder1',filename);
myformat = '%*s';
myformat2 = '%s';
for i = 1 : n_molecules 
    myformat = strcat (myformat, '%f'); 
end
for i = 1 : n_molecules 
    myformat2 = strcat (myformat2, '%*f'); 
end
ar = textscan (myfid, myformat2, -1);
name = ar{1, 1};
name = [name ar{1, 1}];
fclose (myfid);
myfid =  fopen ('result.txt', 'r');
cell_array = textscan (myfid, myformat, -1);
data = cell_array{1, 1};
for i = 2 : n_molecules
    data = [data cell_array{1, i}];
end
fclose (myfid);
n_molecules = n_molek;
%data=data';
size(data)
i1=1;
for i=1:MAX
    if(i1>MAX)
        break;
    end
    j1=i1+1;
    for j=i+1:MAX
        if(j1>MAX)
            break;
        end
        cor=coreluacua (data(i1,:),data(j1,:));
        if cor==0
            data=[data(1:j1-1,:);data(j1+1:MAX,:)];
            name=[name(1:j1-1,:);name(j1+1:MAX,:)];
            j1=j1-1;
            MAX=MAX-1;
        end
        j1=j1+1;
    end
    i1=i1+1;
end
%data=[data(1:168-1,:);data(168+1:MAX,:)];
%name=[name(1:168-1,:);name(168+1:MAX,:)];
%MAX=MAX-1;
size(data)
mkdir ('out');
fid = fopen (fullfile ('out','result1.txt'), 'w') ;
for j=1:MAX
   name1= char(name(j));
   fprintf (fid,'%s               ',name1);
   for i=1:n_molecules
        fprintf (fid,'%f ',data(j,i));
   end
   fprintf (fid,'\r\n');
end
fclose(fid);
%data = textread(fullfile (matlabroot, 'work', 'ћ√”јk-NN',filename));
y = textread('y.txt');
y1=y(:,1);
y=y1;
size(y)
sampl = data;
for i=1:n_molecules
    for j=1:kol_coced
        Class(i,j)=0;
    end
end
for i=1:5
    for j=1:kol
        result(i,j)=0;
    end
end
num=0;
max=0;
J=0;
k1=0;
t=1;
j1=0;
filename='result.txt';
fid = fopen (fullfile ('out',filename), 'a') ;
for j=1:MAX
    
    sampl1=sampl(j,:);
    if Rad~=0
        raduys=opr_raduys(sampl1,Rad);
    else
        raduys=10000000000;
    end
    B=sampl1(2:n_molecules);
    y2=y(2:n_molecules,:);
    training = B';
    group = y2;
    Molekula=sampl1(1,1);
    Class(1,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    for i=2:n_molecules-1
        A=sampl1(1:i-1);
        y1=y(1:i-1,:);    
        B=sampl1(i+1:n_molecules);
        y2=y(i+1:n_molecules,:);
        training = [A,B]';
        group = [y1;y2];
        Molekula=sampl1(1,i);
        if i==33
        end
        Class(i,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    end
    A=sampl1(1:n_molecules-1);
    y1=y(1:n_molecules-1,:);    
    training = A';
    group = y1;
    Molekula=sampl1(1,n_molecules);
    Class(n_molecules,:) = Knnv4_1(Molekula, training, group,kol_coced,M,raduys);
    Class1=Class;
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
        if(n_molecules-otkaz(i)~=0)
            y_sr=y_sr/(n_molecules-otkaz(i));
        end
        for k=1:n_molecules
            if(Class(k,i)==1000000)
                continue;
            end
            R_2(i)= R_2(i)+ (Class(k,i)-y(k))*(Class(k,i)-y(k));
            y_pr =y_pr+(y_sr-y(k))*(y_sr-y(k));
        end
        if (y_pr~=0)
            R_2(i)= 1-(R_2(i)/y_pr);
        else
            R_2(i)=0;
        end
    end
    R_2
   max=-1000000000000000;
    for i1=1:kol_coced
        
        neopr=otkaz(i1);
        Da=0;
   %     for i=1:n_molecules
    %        if Class(i,i1)>l_granutcua&&Class(i,i1)<r_granutcua&&(Class(i,i1)+y(i)>0||Class(i,i1)+y(i)<0)
        %        Da=Da+1;
         %   end
          %  if Class(i,i1)+y(i)==0
           %     neopr=neopr+1;
           % end
       % end
        %if koef==2
        %    count=Da/(n_molecules-neopr)*100;
        %else
       % count=Da+koef*neopr;
       % end
        count=R_2(i1);
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
        result(5,num)=J;        
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
            result(5,I)=J;
        end         
    end
    max=-1000000000000000;
    
end
fclose(fid);
for i=1:kol
    for i1=1:kol-1
        if result(1,i1)< result(1,i1+1)
            s=result(:,i1);
            result(:,i1)= result(:,i1+1);
            result(:,i1+1)=s;
        end
    end
end
filename='out.txt';
fid = fopen (fullfile ('out',filename), 'a') ;
fprintf (fid,'Parametru test8v4:\r\n metruk %s kol_disk %d kol_mol %d kol_otdor %d koef %f l_granutcua %f r_granutcua %f glybuna %d kol_coced %d raduys %d\r\n\r\n',M,MAX,n_molek,kol,koef,l_granutcua,r_granutcua,glybuna,kol_coced,Rad);
for i=1:num
    klast=n_molecules-result(3,i);
    if klast==0
        klast=-1;
    else
        klast=result(2,i)/(n_molecules-result(3,i))*100;
    end
    fprintf (fid,'count %f  max %d  neopr %d    k %d    num1 %d    verno %f\r\n',result(1,i),result(2,i),result(3,i),result(4,i),result(5,i),klast);
end
fprintf (fid,'\r\n\r\n\r\n');
fclose(fid);

max=0;
for j1=1:glybuna
    [result1]=estumate_v4(M,MAX,n_molecules,num,koef,l_granutcua,r_granutcua,glybuna,kol_coced,data,result(:,1),y,Rad);
    for j=2:num
        j
        [res]=estumate_v4(M,MAX,n_molecules,num,koef,l_granutcua,r_granutcua,glybuna,kol_coced,data,result(:,j),y,Rad);
        result1=[result1,res];
    end
    kou=0;
    i1=1;
    for i=1:kol*kol
        if result1(1,i)<= result(1,i1);
            result1(1,i)=0;
        end
        kou=kou+1;
        if(kou/kol>0.9999999999)
            i1=i1+1;
            kou=0;
        end
    end
    for i=1:kol*kol
        for i1=1:kol*kol-1
            if result1(1,i1)< result1(1,i1+1)
                s=result1(:,i1);
                result1(:,i1)= result1(:,i1+1);
                result1(:,i1+1)=s;
            end
        end
    end
    fid = fopen (fullfile ('out',filename), 'a') ;
    A=size(result1);
    N=A(1)-5;
    for i=1:num*num
        fprintf (fid,'count %f  max %d  neopr %d    k %d    num1 %d   ',result1(1,i),result1(2,i),result1(3,i),result1(4,i),result1(5,i));
        for i1=1:N
            fprintf (fid,'num%d %d   ',1+i1,result1(5+i1,i));    
        end
        klast=n_molecules-result1(3,i);
        if klast==0
            klast=-1;
        else
            klast=result1(2,i)/(n_molecules-result1(3,i))*100;
        end
        fprintf (fid,'verno %f\r\n',klast);
    end
    fprintf (fid,'\r\n\r\n\r\n');
    fclose(fid);
    result=result1(:,1:num);
end