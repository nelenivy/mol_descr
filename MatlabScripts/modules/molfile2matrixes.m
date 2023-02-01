function [Coordinates, Bound, atomName] = molfile2matrixes(x)
%returns BoundBlock and AtomBlock
count=1;
a=[];
Coordinates=[];
CoordinatesLine=[];
Bound=[];
BoundLine=[];
fid=fopen(x, 'r');

while count<5 %propuskaem stroki-zagolovki
    line=fgets(fid);
    count=count+1;
end

atomCount=str2num(strcat(line(1), line(2), line(3))); %kolichestvo atomov
boundCount=str2num(strcat(line(4), line(5), line(6))); %kolichestvo svyazey

% --------------- AtomBlock -----------------------------------
atomName=[];
count=1;
for n=1:atomCount
    line=fgets(fid); %schitat sled. stroku
    
    for m=1:3   %coordinares
        
        while isspace(line(count))==1
            count=count+1;
        end
        
        while isspace(line(count))==0
            a=cat(2, a, line(count));
            count=count+1;
        end
    
        CoordinatesLine=cat(2, CoordinatesLine, str2double(a));
        a=[];
    
        while isspace(line(count))==1
              count=count+1;
        end
        
        
    end
    name=[];
    
    name=cat(2, name, line(count));
    name=cat(2, name, line(count+1));
    
    atomName=cat(1, atomName, name);
    
    Coordinates=cat(1, Coordinates, CoordinatesLine);
    CoordinatesLine=[];
    count=1;
end

% ----------------------- BondBlock ----------------------------
for n=1:boundCount
    line=fgets(fid); %schitat sled. stroku
    
    for m=1:3
        a=str2num(strtrim(strcat(line(1:3))));
        line=line(4:length(line));
        BoundLine=cat(2, BoundLine, a);
    end
    Bound=cat(1, Bound, BoundLine);
    BoundLine=[];
    count=1;
end
fclose(fid);