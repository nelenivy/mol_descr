function descMatrix=Descriptor(FileName, k, format)
%returns list of coded chains

descMatrix=struct('chain','', 'count', 0);
s=0;

ind=find(FileName=='\',1,'last');
if (fopen(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\',format,'\chains_coded.mat','r'))==-1)
    ChainCode(FileName, k ,format);
    fclose('all');
end
load(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\',format,'\chains_coded.mat'));

if (isempty(chains_coded))
    ind=find(FileName=='\',1,'last');
    warning off
    mkdir(strcat(FileName(1:ind-1),'\linear_fragments'));
    warning on
    NewFileName=strcat(FileName(1:ind-1), '\linear_fragments\level1\',num2str(k),'\',format,'\', FileName(ind+1:find(FileName=='.'-1)));
   
    NewFileName=strcat(NewFileName, 'cc', int2str(k));

    fid=fopen(NewFileName, 'w');
    fclose(fid);
    return;
end

n=size(chains_coded);
descMatrix(s+1).chain=chains_coded(1,:);
descMatrix(s+1).count=1;

for I=2:n(1)
    s=size(descMatrix);
  
    for J=1:s(2)
        ind=findstr(descMatrix(J).chain, chains_coded(I,:));
        if (ind==1)
            NewCount=descMatrix(J).count;
            NewCount=NewCount+1;
            descMatrix(J).count=NewCount;
            break;
        end
    end
        if (isempty(ind))
            descMatrix(s(2)+1).chain=chains_coded(I, :);
            descMatrix(s(2)+1).count=1;
        end
end

descriptors=descMatrix;
clear descMatrix
%----------------- to file -----------------------
ind=find(FileName=='\',1,'last');
save(strcat(FileName(1:ind),'linear_fragments\level1\',num2str(k),'\',format,'\descriptors.mat'),'descriptors');
