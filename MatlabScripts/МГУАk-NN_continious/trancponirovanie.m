Matrux = textread('result.txt');
s=size(Matrux)
fid = fopen ('result.txt', 'w') ;
for j=1:s(2)
    fprintf (fid,'%d ',j);
    for i=1:s(1)
        fprintf (fid,'%d ',Matrux(i,j));
    end
    fprintf (fid,'\r\n');
end
fclose(fid);