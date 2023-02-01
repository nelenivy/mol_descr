function Distance=OT2(dirName, MaxMinMid, Keys)
% считывает не кодированные цепочки из .OT1 
% считывает кодированные цепочки из .codeOT1
% формирует особые точки второго уровня
% записывает в соответствующие файлы
OTMatrix={};
OTCodedMatrix={};
indexStr=0;

Distance=[];
load(strcat(Keys.root_package,'\molecules\',dirName,'\linear_fragments\level1\',...
    Keys.chain_length,'\chains_numb.mat'));         %считывает не кодированные цепочки

load(strcat(Keys.root_package,'\molecules\',dirName,'\linear_fragments\level1\',...
    Keys.chain_length,'\',Keys.markers,'\chains_coded.mat'))  %считывает кодированные цепочки

ChainLength=size(chains_numb,2);
zerstr=zeros(1, ChainLength^2);
DistanseBetweenAtoms=[]; % расстояния между атомами в одной ОТ второго уровня

DistanceMatrix=GeometricDistanceMatrix(strcat(Keys.root_package,'\molecules\',dirName,'\',dirName,'.mol'));% матрица геометрических расстояний между атомами
for I=1:size(chains_numb,1)
    for J=I+1:size(chains_numb,1)
        OTMatrix=cat(1,OTMatrix, [chains_numb(I,:), chains_numb(J,:)]);% матрица не кодированных цепочек
        OTCodedMatrix=cat(1,OTCodedMatrix, strcat(chains_coded(I,:), chains_coded(J,:))); % матрица кодированных цепочек
        DistanseBetweenAtoms=[DistanseBetweenAtoms; zerstr];
        indexCol=1;
        indexStr=indexStr+1;
        for k=1:size(chains_numb,2)
            for l=1:size(chains_numb,2)
                % каждой строке матрицы расстояний соответствует одна пара
                % цепочек
                DistanceBetweenAtoms(indexStr, indexCol)=DistanceMatrix(chains_numb(I,k), chains_numb(J,l)); % заполнение матрицы расстояний
                indexCol=indexCol+1;
            end
        end
    end
end
indexCol=indexCol-1;

% выбор заданного типа расстояния
if(MaxMinMid=='max')
    
    for I=1:indexStr
        maxV=max(DistanceBetweenAtoms(I,:));
        
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, ' ');
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, num2str(maxV));
        Distance=[Distance; maxV];
    end
end

if(MaxMinMid=='min')
    
    for I=1:indexStr
        minV=min(DistanceBetweenAtoms(I,:));
        
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, ' ');
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, num2str(minV));
        Distance=[Distance; minV];
    end
end

if(MaxMinMid=='mid')
    
    for I=1:indexStr
        mid=DistanceBetweenAtoms(I,1);
        for J=1:indexCol
            mid=mid+DistanceBetweenAtoms(I, J);
        end
        mid=mid/indexCol;
        
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, ' ');
        OTCodedMatrix{I}=cat(2, OTCodedMatrix{I}, num2str(mid));
        Distance=[Distance; mid];
    end
end

% записать в файлы
clear chains_number
chains_number=[];
for I=1:size(OTMatrix,1)
    chains_number=[chains_number; OTMatrix{I}];
end
clear OTMatrix;
warning off
mkdir(strcat(Keys.root_package,'\molecules\',dirName,'\',Keys.method,'\', Keys.chain_length,'\',Keys.markers));
warning on
save(strcat(Keys.root_package,'\molecules\',dirName,'\',Keys.method,'\', Keys.chain_length,'\chains_number.mat'),'chains_number');


%сделать общую матрицу, убрать повторения
clear chains_coded
clear descriptors
chains_distance=sort(OTCodedMatrix);


save(strcat(Keys.root_package,'\molecules\',dirName,'\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,'\chains_distance.mat'),'chains_distance');

