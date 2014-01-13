function [] = claster(num_dick,kol_kl,n_molek)
n_molecules=n_molek;
kol_kl=kol_kl-1;
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
y = textread(fullfile (matlabroot, 'work', 'МГУАk-NN','y.txt'));
y1=y(:,1);
y=y1;
s=size(data)
if num_dick==0
    sz=s(1);
else
    sz=0;
end
for i=1:n_molek
    num_mol(i)=i;
end
e=0.0000000001;
for i=1:s(2)
    min(i)=10000000;
    J(i)=0;
end

for i=1:s(2)
    for j=1:s(2)
        x=0;
        if(j~=i)
            for k=1:sz
                x=x+(data(k,i)-data(k,j))*(data(k,i)-data(k,j));
            end
            if num_dick~=0
                x=x+(data(num_dick,i)-data(num_dick,j))*(data(num_dick,i)-data(num_dick,j));
            end  
            if(min(i)>x)
                min(i)=x;
                J(i)=j;
            end
        end
    end
end
I=1;
MIN=min(1);
for i=2:s(2)
    if(min(i)<MIN)
        MIN=min(i);
        I=i;
    end
end
derevo(1,1)=I;
derevo(2,1)=0;
derevo(1,2)=J(I);
derevo(2,2)=MIN;
for i=1:s(1)
    derevo(2+i,1)=data(i,I);
    derevo(2+i,2)=data(i,J(I));
end
derevo(3+s(1),1)=num_mol(I);
derevo(3+s(1),2)=num_mol(J(I));
derevo(4+s(1),1)=y(I);
derevo(4+s(1),2)=y(J(I));
if(derevo(1,1)<derevo(1,2))
    data=[data(:,1:derevo(1,1)-1),data(:,derevo(1,1)+1:derevo(1,2)-1),data(:,derevo(1,2)+1:s(2))];
    num_mol=[num_mol(1:derevo(1,1)-1),num_mol(derevo(1,1)+1:derevo(1,2)-1),num_mol(derevo(1,2)+1:s(2))];
    y=[y(1:derevo(1,1)-1);y(derevo(1,1)+1:derevo(1,2)-1);y(derevo(1,2)+1:s(2))];
else
    data=[data(:,1:derevo(1,2)-1),data(:,derevo(1,2)+1:derevo(1,1)-1),data(:,derevo(1,1)+1:s(2))];
    num_mol=[num_mol(1:derevo(1,2)-1),num_mol(derevo(1,2)+1:derevo(1,1)-1),num_mol(derevo(1,1)+1:s(2))];
    y=[y(1:derevo(1,2)-1);y(derevo(1,2)+1:derevo(1,1)-1);y(derevo(1,1)+1:s(2))];
end
n=s(2)-2;
for i1=1:n
    s=size(data);
    d=size(derevo);
    for i=1:s(2)
        min(i)=10000000000;
    end
    for i=1:s(2)
        for j=1:d(2)
            X(j)=0;
                for k=1:sz
                    X(j)=X(j)+(data(k,i)-derevo(2+k,j))*(data(k,i)-derevo(2+k,j));
                end
                if num_dick~=0
                    X(j)=X(j)+(data(num_dick,i)-derevo(2+num_dick,j))*(data(num_dick,i)-derevo(2+num_dick,j));
                end
                if(min(i)>X(j))
                    min(i)=X(j);
                end
        end
    end
    I=1;
    MIN=min(1);
    for i=2:s(2)
        if(min(i)<MIN)
            MIN=min(i);
            I=i;
        end
    end
    derevo(1,i1+2)=I;
    derevo(2,i1+2)=MIN;
    for i=1:s(1)
        derevo(2+i,i1+2)=data(i,I);
    end
    derevo(3+s(1),i1+2)=num_mol(I);
    derevo(4+s(1),i1+2)=y(I);
    data=[data(:,1:derevo(1,i1+2)-1),data(:,derevo(1,i1+2)+1:s(2))];
    num_mol=[num_mol(1:derevo(1,i1+2)-1),num_mol(derevo(1,i1+2)+1:s(2))];
    y=[y(1:derevo(1,i1+2)-1);y(derevo(1,i1+2)+1:s(2))];
end
size(derevo) 
for i=1:kol_kl
    max=0;
    for j=1:n+2
        if derevo(2,j)>max
            max=derevo(2,j);
        end
    end
    for j=1:n+2
        if (derevo(2,j)- max < e)&&(derevo(2,j)- max > -e)
            derevo(2,j)=-1;
        end
    end
end
num=1;
mkdir (fullfile (matlabroot, 'work','МГУАk-NN','out'));
cd('out');
OutFile='claster';
OutFile=cat(2, OutFile, num2str(num));
OutFile=cat(2, OutFile, '.txt');
OutId=fopen(OutFile, 'w');
OutFile1='claster_inf';
OutFile1=cat(2, OutFile1, '.txt');
OutId1=fopen(OutFile1, 'w');
OutFile2='y';
OutFile2=cat(2, OutFile2, num2str(num));
OutFile2=cat(2, OutFile2, '.txt');
OutId2=fopen(OutFile2, 'w');
for i=1:n+2
    if derevo(2,i)+1==0
        fprintf (OutId1,'\r\n');
        fclose(OutId);
        fclose(OutId2);
        num=num+1;
        OutFile='claster';
        OutFile=cat(2, OutFile, num2str(num));
        OutFile=cat(2, OutFile, '.txt');
        OutId=fopen(OutFile, 'w');
        OutFile2='y';
        OutFile2=cat(2, OutFile2, num2str(num));
        OutFile2=cat(2, OutFile2, '.txt');
        OutId2=fopen(OutFile2, 'w');
        for j=1:s(1)
            fprintf (OutId,'%d   ',derevo(j+2,i));
        end
        fprintf (OutId1,'%d   ',derevo(s(1)+3,i));
        fprintf (OutId2,'%d   ',derevo(s(1)+4,i));
        fprintf (OutId2,'\r\n');
        fprintf (OutId,'\r\n');
    else
        for j=1:s(1)
            fprintf (OutId,'%d   ',derevo(j+2,i));
        end
        fprintf (OutId1,'%d   ',derevo(s(1)+3,i));
        fprintf (OutId2,'%d   ',derevo(s(1)+4,i));
        fprintf (OutId2,'\r\n');
        fprintf (OutId,'\r\n');
    end
end
fprintf (OutId1,'\r\n');
fprintf (OutId1,'%d   условие кластеризации',sqrt(max)/3);
fprintf (OutId1,'\r\n');
fclose(OutId1);
fclose(OutId);
fclose(OutId2);
cd('..');
           
