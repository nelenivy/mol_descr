function cluster_MST(KeyFile, YFile)

% функция кластеризации матрицы методом минимального покрывающего дерева. Результат сохраняется в
% директорию с полным путем, соответствующем KeyFile
% minimum spanning tree

Keys=readKeys(KeyFile);
path=strcat(Keys.method,'\',Keys.chain_length,'\',Keys.markers);
fid=fopen(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'),'r');
if(fid<0)
    error('MDMatrix is not exist');
end
fclose(fid);
load(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'));
y=load(YFile);