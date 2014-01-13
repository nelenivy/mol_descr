function [res] = raduys_grafa(data)
s=size(data);
for i=1:s(2)
    MAX(i)=0;
end
for i=1:s(2)
    for j=1:s(2)
        x=0;
        for k=1:s(1)
            x=x+(data(k,i)-data(k,j))*(data(k,i)-data(k,j));
        end
        if(MAX(i)<x)
            MAX(i)=x;
        end 
    end
end
for i=1:s(2)-1
    if(MAX(i)<MAX(i+1))
        MAX(i+1)=MAX(i);
    end
end
res=MAX(s(2));

            
