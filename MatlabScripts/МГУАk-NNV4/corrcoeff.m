function c=corrcoeff(y1,y2,dim)
c=0;

for i=1:dim
    c=c+(y1(i)-mean(y1))*(y2(i)-mean(y2));
end

c1=0;
for i=1:dim
    c1=c1+(y1(i)-mean(y1))*(y1(i)-mean(y1));
end

c2=0;
for i=1:dim
    c2=c2+(y2(i)-mean(y2))*(y2(i)-mean(y2));
end

c=c/sqrt(c1*c2);