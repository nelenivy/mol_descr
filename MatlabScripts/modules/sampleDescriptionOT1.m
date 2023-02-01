function sampleDescriptionOT1(KeyFile1, KeyFile2, SDFFile)

% метод для построения МД-матриц для новой выборки молекул с неизвестной активностью
% На вход подается 2 файла с ключами. 
% KeyFile1 - параметры хранения результатов, 
% по KeyFile2 определяются списки дескрипторов, с помощью которых будут
% строиться МД-матрицы

display('sampleDescriptionOT1...')
warning off;
if nargin<3
    SDFFile='0';
end
Keys1=readKeys(KeyFile1);
    Keys2=readKeys(KeyFile2);
if (~strcmp(SDFFile,'0'))
    cutSDF(SDFFile, KeyFile1);
else
    
    if (~isdir(strcat(Keys1.root_package)))
    end
end
connectivity_profile(KeyFile1); 
Descriptor_all(KeyFile1);
Descriptor_list(KeyFile1);
 fid=fopen(strcat(pwd,'\',Keys2.root_package,'\descriptor_list\',Keys1.method,'\', Keys1.chain_length,'\',Keys1.markers,'\descriptor_list.mat'),'r');
 if(fid>0)
     load(strcat(pwd,'\',Keys2.root_package,'\descriptor_list\',Keys1.method,'\', Keys1.chain_length,'\',Keys1.markers,'\descriptor_list.mat'));
     MDmatrix(KeyFile1, descriptor_list);
 else
     error('description is not exist \n');
 end