function Distance=OT2_all(KeyFile, DistFile)
% построение особых точек второго уровен€ linear_fragments\level 2 дл€ всех
% молекул в root_package\molecules и префиксом mol_prefix
% ¬ход: файл с параметрами запуска KeyFile, файл с заданным рассто€нием и интервалами 
% ¬ыход: вектор всех попарных рассто€ний между цепочками молекулы

display('OT2_all...');
Keys=readKeys(KeyFile);
if((isempty(Keys.root_package)) || (isempty(Keys.mol_prefix)) ||(isempty(Keys.chain_length))|| (isempty(Keys.markers)) || (isempty(Keys.method)))
    return;
end
% если на текущий момент не построены ќ“ 1 уровн€ с заданными параметрами
if (fopen(strcat(Keys.root_package,'\molecules\',Keys.mol_prefix,'0001\linear_fragments\level1\',Keys.chain_length,'\',Keys.markers,'\chains_coded.mat'),'r')<0)
    %error('Build the descriptors of 1st level');
    fid=fopen('tmp.txt','w');
    field=fields(Keys);
    for I=1:size(field,1)
        if (~strcmp(field{I},'method'))
            fprintf(fid, strcat(field{I},'=',Keys.(field{I}),'\n'));
        else
            fprintf(fid, strcat(field{I},'=linear_fragments\\level1\n'));
        end
    end
    fclose(fid);
    OT1main('tmp.txt');
    delete('tmp.txt');
end

choice=1;
if (fopen(strcat(Keys.root_package,'\molecules\',Keys.mol_prefix,'0001\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,'\chains_distance.mat'),'r')>0)
    choice=input('OT2_all.m The chains are already builded. Press the key \n 1. Rebuild \n 2. Use the builded chains\n');
end
load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));
[MaxMinMid,intervals]=readIntervals(DistFile);
Distance=[];
files=dir(strcat(pwd,'\', Keys.root_package,'\molecules')); %imena fseh direktoriy
n=size(files);
for I=3:n(1)
    k=findstr(files(I).name, Keys.mol_prefix);
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1) && isempty(find(strcmp(exception,files(I).name))))
            NewDir=strcat(pwd,'\',Keys.root_package,'\molecules\',files(I).name);
            fileName=strcat(NewDir,'\',files(I).name, '.mol'); % formiruem im9 nujnogo fayla
            fid=fopen(fileName,'r');
            if (fid>0 && choice==1)
               Distance=[Distance; OT2(files(I).name,MaxMinMid,Keys)];
            end
            fclose(fid);
        end
    end
end
