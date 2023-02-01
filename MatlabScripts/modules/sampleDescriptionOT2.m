function sampleDescriptionOT2(KeyFile1, KeyFile2, DistFile, SDFFile)

% метод для построения МД-матриц для новой выборки молекул с неизвестной активностью
% На вход подается 2 файла с ключами. 
% KeyFile1 - параметры хранения результатов, 
% по KeyFile2 определяются списки дескрипторов, с помощью которых будут
% строиться МД-матрицы

display('sampleDescriptionOT2...')
warning off;
if nargin<4
    SDFFile=0;
end
Keys1=readKeys(KeyFile1);
    Keys2=readKeys(KeyFile2);

if (SDFFile~=0)
    cutSDF(SDFFile, KeyFile1);
end
Distance=OT2_all(KeyFile1, DistFile);
out=CodeAllOT(KeyFile1, DistFile);
fid=fopen(strcat(pwd,'\',Keys2.root_package,'\descriptor_list\',Keys1.method,'\', Keys1.chain_length,'\',Keys1.markers,'\descriptor_list.mat'),'r');
if(fid>0)
    load(strcat(pwd,'\',Keys2.root_package,'\descriptor_list\',Keys1.method,'\', Keys1.chain_length,'\',Keys1.markers,'\descriptor_list.mat'));
    if (out==0)
        Descriptor_list(KeyFile1);
        MDmatrix(KeyFile1,descriptor_list);
    end
else
    error('description is not exist \n');
end