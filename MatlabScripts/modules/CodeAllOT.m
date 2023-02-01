function output=CodeAllOT(KeyFile, DistFile)
display('CodeAllOT...');
% замена численного значения расстояния между цепочками в молекуле на имена
% интервалов. Интервалы задаются в DistFile
output=0;
[MaxMinMid,intervals]=readIntervals(DistFile);
Keys=readKeys(KeyFile);
load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));

files=dir(strcat(Keys.root_package,'\molecules'));
if (isempty(intervals))
    warning('intervals are not defined');
    output=1;
    return;
end
for I=3:size(files,1)
    k=findstr(files(I).name, Keys.mol_prefix);
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1) && isempty(find(strcmp(exception,files(I).name))))
            load(strcat(Keys.root_package,'\molecules\',files(I).name,'\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,...
                '\chains_distance.mat'))
            
            for J=1:size(chains_distance,1)
                distance=str2num(strtrim(chains_distance{J}(5*str2num(Keys.chain_length)*2+1:length(chains_distance{J}))));
                for l=1:size(intervals,1)
                    if ((distance>=intervals{l,2}) && (distance<=intervals{l,3}))
                        distance=intervals{l,1};
                        break;
                    end
                end
                chains_distance{J}=cat(2,chains_distance{J}(1:5*str2num(Keys.chain_length)*2+1),distance);
            end
            chains_distance=sort(chains_distance);
            clear chains_coded;
            clear descriptors;
            descriptors.chain=chains_distance{1};
            descriptors.count=1;
            chains_coded{1}=chains_distance{1};
            for J=2:size(chains_distance,1)
                ind=binary_search(chains_coded, chains_distance{J});
                if (ind==0)
                    chains_coded=cat(1,chains_coded,chains_distance{J});
                    descriptors(size(descriptors,2)+1).chain=chains_distance{J};
                    descriptors(size(descriptors,2)).count=1;
                else
                    descriptors(ind).count=descriptors(ind).count+1;
                end 
            end
           
            chainMatrix=[];
            for J=1:size(chains_coded,1)
                chainMatrix=[chainMatrix;chains_coded{J}];
            end
            chains_coded=chainMatrix;
            clear chainMatrix;
            save(strcat(Keys.root_package,'\molecules\',files(I).name,'\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,...
                '\chains_coded.mat'),'chains_coded');
            save(strcat(Keys.root_package,'\molecules\',files(I).name,'\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,...
                '\descriptors.mat'),'descriptors');
        end
    end
end

