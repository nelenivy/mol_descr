function res=TDistance(ConnectivMatrix, StartInd, FinishInd)
% topologi4eskoe rasstoyanie mejdu StartInd i FinishInd v ConnectivMatrix
T=0;
OldFront=[];
NewFront=[];
Mark=[];
count=size(ConnectivMatrix);
for I=1:count(1)
    Mark(I)=-1;
end
Mark(StartInd)=0;
OldFront=cat(2, OldFront, StartInd);

OFsize=size(OldFront);
while OFsize(2)~=0
		for J=1:OFsize(2)		%proyti po matrice svyaznosti
			for I=1:count(1)
				if (I==OldFront(J))
					continue;
                end

				%nashli incid. neotmechen. vershinu
				if (ConnectivMatrix(OldFront(J), I)~=0) && (Mark(I)==-1)
					NewFront=cat(2, NewFront, I);
					Mark(I)=T+1;
                end
            end
			
			%esli vstretili nujnuyu vershinu, vernut rezultat
            NFsize=size(NewFront);
			for I=1:NFsize(2)
				if (FinishInd==NewFront(I))
					res=T+1;
                    return;
                end
            end

        end
		%esli ne nashli vershinu, sled. iteraciya
		T=T+1;
		OldFront=[];
		for I=1:NFsize(2)
			OldFront=cat(2, OldFront, NewFront(I));
        end
        
		NewFront=[];
OFsize=size(OldFront);

end

res=0;