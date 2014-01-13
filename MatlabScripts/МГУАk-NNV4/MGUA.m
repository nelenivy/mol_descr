function R=MGUA(mol,desk,file_md,file_activity,glubina,q,Rt)
f=fopen(file_md,'r');
%for i=1:mol
    %for j=1:desk
        MD_matrix=fscanf(f,'%f',[mol desk]);
    %end
%end
fclose(f);

f=fopen(file_activity,'r');
%for i=1:mol
    %for j=1:desk
        activity=fscanf(f,'%f',[mol 1]);
        size(activity);
    %end
%end
fclose(f);

for i=1:desk
    a=mean(MD_matrix(1:mol,i));
    for j=1:mol
        MD_matrix(j,i)=MD_matrix(j,i)-a;
    end
end
MD_matrix;
a=mean(activity);

for j=1:mol
    activity(j)=activity(j)-a;
end
activity
%MD_matrix=textread(file,'%d')
[m,n]=size(MD_matrix);
count=1;
%����� ������ ��� ��������
for i=1:n
    for j=i:n        
        [x(count,1:2),x(count,3)]=lsqnonneg([MD_matrix(1:m,i),MD_matrix(1:m,j)],activity);
        x(count,4)=i;
        x(count,5)=j;
        count=count+1;
    end
end

max=x(1,1:5);
maxj=1;
y(1:q,1:5)=0;
%�������������� �� ������
for i=1:q
    for j=1:(count-1)
        if x(j,3)<max(1,3)
            max(1,1:5)=x(j,1:5);
            maxj=j;
        end
    end
    y(i,1:5)=max(1,1:5);
    x(maxj,3)=max(3)*2;
end

err1=y(1,3);

for i=1:q
    for j=1:m
    Q(j,i)=y(i,1)*MD_matrix(j,y(i,4))+y(i,2)*MD_matrix(j,y(i,5));
    end
end


err2=err1;
%���� ������ ���������
for t=1:glubina
    count=1;
    
    for i=1:q
        for j=1:n
            [z(count,1:2),z(count,3)]=lsqnonneg([Q(1:m,i),MD_matrix(1:m,j)],activity);
            z(count,4)=i;
            z(count,5)=j;
            count=count+1;
        end
    end
    count=count-1;    
    y(1:count,1:5)=0;

    for i=1:count
        max(1,1:5)=z(1,1:5);
        maxj=1;
        for j=1:(count)
            if z(j,3)<max(1,3)
                max(1,1:5)=z(j,1:5);
                maxj=j;
            end
        end
        y(i,1:5)=max(1,1:5);
        z(maxj,3)=2000000000*max(1,3);
    end
    
    for i=1:m
        Q1(i,1)=y(1,1)*Q(i,y(1,4))+y(1,2)*MD_matrix(i,y(1,5));
    end    
    
    j=2;
    i=2;
        while (i<=q)&&(j<=count)
            corr=1;
            j;
            for k=1:m
                c(k)=y(j,1)*Q(k,y(j,4))+y(j,2)*MD_matrix(k,y(j,5));
            end
            
            for k=1:i-1
                corrcoeff(c,Q1(1:m,k),m);
                if corrcoeff(c,Q1(1:m,k),m)>Rt
                    corr=0;
                end
            end
            
            if corr==1
               j ;
               for k=1:m
                    Q1(k,i)=c(k);
               end 
               
               i=i+1;
            end
            
            j=j+1;
        end
    
    err1=err2
    err2=y(1,3)
   
    for i=1:m
        for j=1:q
            Q(i,j)=Q1(i,j);
        end
    end
    
end

y1=0;
%������� ��-��������
for i=1:m
  y1=y1+(activity(i)-mean(activity))*(activity(i)-mean(activity));
end
R=1-err2*err2/y1;

f=fopen('result.txt','w');
fprintf(f,'R=%f err=%f',R,err2);
fclose(f);