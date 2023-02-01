function chainsNew=ReversEqual(matrix)
%delete repeted chains like a-b-c; c-b-a from matrix
n=size(matrix);
chainsNew=[];

for I=1:n(1)
    flag=0;
    for J=I+1:n(1)
        for k=1:n(2)
            if (matrix(I,k)==matrix(J,n(2)-k+1))
                flag=1;
                                    
            else
                flag=0;
                break;
            end
        end
        
        if (flag==1)
            break;
        end
    end
    
    if (flag==0)
        chainsNew=cat(1, chainsNew, matrix(I,:));
    end
end
