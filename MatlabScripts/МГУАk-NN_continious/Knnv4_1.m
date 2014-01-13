function [Res] = Knnv4_1(Molekula, training, group,kol,M,raduys);
e=0.0000000001;
A=size(Molekula);
siz1=A(2);
A=size(training);
siz2=A(1);
sosed=0;
for k=1 : kol
    dist(k)=0;
    pres(k)= 0;
    Res(k)= 0;
    count(k)=0;
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
                 pres(l)= pres(l)+group(j);
                 count(l)= count(l)+1;
                 f=1;
               end
           end
           if f==0
                sosed=sosed+1;
                dist(sosed)=s;
                pres(sosed)=group(j);
                count(sosed)=1;
           end
        else
           f=0;
           for l=1 :sosed
               if s - dist(l)<e && s - dist(l)>-e
                 count(l)= count(l)+1;
                 pres(l)= pres(l)+group(j);
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
                if s - dist(num)<e && s - dist(num)>-e
                      pres(num)= pres(num)+group(j);
                      count(num)= count(num)+1;
                end
                if s - dist(num)<-e
                    dist(num)=s;
                    pres(num)=group(j);
                    count(num)=1;
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
                s=pres(i);
                pres(i)=pres(i+1);
                pres(i+1)=s;
                s=count(i);
                count(i)=count(i+1);
                count(i+1)=s;
            end
        end
    end
for k=1 : kol
    if k==500
        fyfh=1;
    end
    kol_s=0;
    for j=1 : k
        if(dist(j)<raduys)
            kol_s=kol_s+count(j);    
        end
    end
    for j=1 : k
        if(dist(j)<raduys)
            Res(k)=Res(k)+pres(j);
        end
    end
    if(kol_s==0)
        Res(k)=1000000;
    else
        Res(k)=Res(k)/kol_s;
    end
end

