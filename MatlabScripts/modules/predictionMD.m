function predictionMD(KeyFile1, KeyFile2)

% предсказание активности для новых соединений
% путь к МД-матрицам новых соединений находится из KeyFile1
% результаты МГУА в KeyFile2

Keys1=readKeys(KeyFile1);
Keys2=readKeys(KeyFile2);

load(strcat(Keys1.root_package,'\MDmatrix\',Keys1.method,'\',Keys1.chain_length,'\',Keys1.markers,'\matrix.mat'));
load(strcat(Keys2.root_package,'\MGUA\MDmatrix\',Keys1.method,'\',Keys1.chain_length,'\',Keys1.markers,'\MGUA.mat'));

% удалить столбцы excluded
ReducedMatrix=zeros(size(matrix,1),1);
if (excluded(1)~=0)
    for I=1:size(matrix,2)
        if (isempty(find(excluded==I)))
            ReducedMatrix=[ReducedMatrix, matrix(:,I)];
        else
            ReducedMatrix(:,1)=ReducedMatrix(:,1) + matrix(:,I);
        end
    end
end
activity=[];
% вычислить активность
for I=1:size(ReducedMatrix,1)
    activity=[activity; sign(coefficients* transpose([ReducedMatrix(I, informative), 1]))];
end

% сохранить значение активности

warning off;
mkdir(strcat(Keys1.root_package,'\predicted_activity\MDmatrix\',Keys1.method,'\',Keys1.chain_length,'\',Keys1.markers));
warning on;
save(strcat(Keys1.root_package,'\predicted_activity\MDmatrix\',Keys1.method,'\',Keys1.chain_length,'\',Keys1.markers,'\activity.mat'),'activity');

