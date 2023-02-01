function [distance, intervals]=readIntervals(DistFile)
% vhodnoy file soderjit tip rasstoyaniaya distance=max/min/mid
% i intervaly v formate a=1:2    \n
distance='';
intervals={};
IntervalsStr=0;

fid=fopen(DistFile, 'r');
if (fid==-1)
    return;
end
line=fgets(fid);
while (line~=-1)
    line=strtrim(line);
    ind=findstr(line,'//');
    if (~isempty(ind))
        NewLine=line(1:ind(1));
    else
        NewLine=line;
    end
    
    if (~isempty(NewLine))
        length=size(NewLine);
        distInd=findstr(NewLine, 'distance');
        if (~isempty(distInd))
            eq=findstr(NewLine, '=');
            if (~isempty(eq))
                distance=NewLine(eq+1:length(2));
                distance=strtrim(distance);
            end
        else
            eq=findstr(NewLine, '=');
            if (~isempty(eq))
                str=NewLine(1:eq-1);
                IntervalsStr=IntervalsStr+1;
                intervals{IntervalsStr,1}=str;
                dp=findstr(NewLine, ':');
                if(~isempty(dp))
                    str=NewLine(eq+1:dp-1);
                    str=strtrim(str);
                    intervals{IntervalsStr,2}=str2num(str);
                    str=NewLine(dp+1:length(2));
                    str=strtrim(str);
                    intervals{IntervalsStr,3}=str2num(str);
                end
            end
        end
    end
    line=fgets(fid);
end
