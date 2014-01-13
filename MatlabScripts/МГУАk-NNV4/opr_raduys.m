function [raduys] = opr_raduys(data,Rad)
s=size(data);
for i=1:s(2)
    min(i)=10000000;
    J(i)=0;
end

for i=1:s(2)
    for j=1:s(2)
        x=0;
        if(j~=i)
            for k=1:s(1)
                x=x+(data(k,i)-data(k,j))*(data(k,i)-data(k,j));
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
if(derevo(1,1)<derevo(1,2))
    data=[data(:,1:derevo(1,1)-1),data(:,derevo(1,1)+1:derevo(1,2)-1),data(:,derevo(1,2)+1:s(2))];
else
    data=[data(:,1:derevo(1,2)-1),data(:,derevo(1,2)+1:derevo(1,1)-1),data(:,derevo(1,1)+1:s(2))];
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
                for k=1:s(1)
                    X(j)=X(j)+(data(k,i)-derevo(2+k,j))*(data(k,i)-derevo(2+k,j));
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
    data=[data(:,1:derevo(1,i1+2)-1),data(:,derevo(1,i1+2)+1:s(2))];
end
%S=0;
%for i=2:n+2
 %   S=S+derevo(2,i);
%end
%if n+1~=0
 %   raduys=S/(n+1);
%else
 %   raduys=0;
%end
for i=1:Rad
    for j=1:n+1
        if derevo(2,j)>derevo(2,j+1)
            w=derevo(2,j+1);
            derevo(2,j+1)=derevo(2,j);
            derevo(2,j)=w;
        end
    end,
end
raduys=derevo(2,Rad);
           
