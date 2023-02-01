function Distance=OT2(dirName, MaxMinMid, Keys)
% ��������� �� ������������ ������� �� .OT1 
% ��������� ������������ ������� �� .codeOT1
% ��������� ������ ����� ������� ������
% ���������� � ��������������� �����
OTMatrix={};
OTCodedMatrix={};
indexStr=0;

Distance=[];
load(strcat(Keys.root_package,'\molecules\',dirName,'\linear_fragments\level1\',...
    Keys.chain_length,'\chains_numb.mat'));         %��������� �� ������������ �������

load(strcat(Keys.root_package,'\molecules\',dirName,'\linear_fragments\level1\',...
    Keys.chain_length,'\',Keys.markers,'\chains_coded.mat'))  %��������� ������������ �������

ChainLength=size(chains_numb,2);
zerstr=zeros(1, ChainLength^2);
DistanseBetweenAtoms=[]; % ���������� ����� ������� � ����� �� ������� ������

DistanceMatrix=GeometricDistanceMatrix(strcat(Keys.root_package,'\molecules\',dirName,'\',dirName,'.mol'));% ������� �������������� ���������� ����� �������
for I=1:size(chains_numb,1)
    for J=I+1:size(chains_numb,1)
        OTMatrix=cat(1,OTMatrix, [chains_numb(I,:), chains_numb(J,:)]);% ������� �� ������������ �������
        OTCodedMatrix=cat(1,OTCodedMatrix, strcat(chains_coded(I,:), chains_coded(J,:))); % ������� ������������ �������
        DistanseBetweenAtoms=[DistanseBetweenAtoms; zerstr];
        indexCol=1;
        indexStr=indexStr+1;
        for k=1:size(chains_numb,2)
            for l=1:size(chains_numb,2)
                % ������ ������ ������� ���������� ������������� ���� ����
                % �������
                DistanceBetweenAtoms(indexStr, indexCol)=DistanceMatrix(chains_numb(I,k), chains_numb(J,l)); % ���������� ������� ����������
                indexCol=indexCol+1;
            end
        end
    end
end
indexCol=indexCol-1;

% ����� ��������� ���� ����������
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

% �������� � �����
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


%������� ����� �������, ������ ����������
clear chains_coded
clear descriptors
chains_distance=sort(OTCodedMatrix);


save(strcat(Keys.root_package,'\molecules\',dirName,'\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,'\chains_distance.mat'),'chains_distance');

