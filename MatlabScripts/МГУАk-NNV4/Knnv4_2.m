function [Res1] = Knnv4_1(Molekula, training, group,kol,M,raduys);
e=0.0000000001;
A=size(Molekula);
siz1=A(2);
A=size(training);
siz2=A(1);
sosed=0;
group_size=4;
for k=1 : kol
    dist(k)=0;
    for j=1:group_size
        pres(k,j)= 0;
        Res(k,j)= 0;
    end 
    Res1(k)=0;
end
for j = 1 : siz2
    s=0;
    for i = 1 : siz1
        s=s+(Molekula(i)-training(j,i))*(Molekula(i)-training(j,i));
    end
        if sosed<kol
           f=0;
           for l=1 :sosed
               if s - dist(l)<e && s - dist(l)>-e
                 pres(l,group(j))= pres(l,group(j))+1;
                 f=1;
               end
           end
           if f==0
                sosed=sosed+1;
                dist(sosed)=s;
                pres(sosed,group(j))=1;
           end
        else
           f=0;
           for l=1 :sosed
               if s - dist(l)<e && s - dist(l)>-e
                 pres(l,group(j))= pres(l,group(j))+1;
                 f=1;
               end
           end
           if f==0
                num=1;
                max=dist(1);
                for i=1:sosed-1
                    if max<dist(i+1)
                        num=i+1;
                        max=dist(i+1);
                    end
                end
                if s - dist(num)<-e
                    dist(num)=s;
                    for i1=1:group_size
                        pres(num,i1)= 0;
                    end
                    pres(num,group(j))=1;
                end
           end
        end        
end
    for j=1 : sosed
        for i=1 : sosed-j
            if dist(i)>dist(i+1)
                s=dist(i);
                dist(i)=dist(i+1);
                dist(i+1)=s;
                s=pres(i,:);
                pres(i,:)=pres(i+1,:);
                pres(i+1,:)=s;
            end
        end
    end
for k=1 : kol
    for j=1 : k
        if(dist(j)<raduys)
            Res(k,:)=Res(k,:)+pres(j,:);
        end
    end
end
for k=1 : kol
    max=0;
    num=0;
    for j=1 : group_size
        if(max<Res(k,j))
            max=Res(k,j);
            num=j;
        end 
    end
    if(num>0.1)
        Res(k,num)=0;
    end
    for j=1 : group_size
        if(max==Res(k,j))
            num=0;
        end 
    end
    Res1(k)=num;
end

