function cluster_kmeans(KeyFile, YFile)
% функция кластеризации матрицы методом к средних. Результат сохраняется в
% директорию с полным путем, соответствующем KeyFile
display('Clustering...');
warning off;
Keys=readKeys(KeyFile);
path=strcat(Keys.method,'\',Keys.chain_length,'\',Keys.markers);
fid=fopen(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'),'r');
if(fid<0)
    error('MDMatrix is not exist');
end
fclose(fid);
load(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'));
y=load(YFile);
% обработка несвязных молекул
load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));
if(~isempty(exception))
    numexc=[];
    for I=1:size(exception,1)
        numexc=[numexc; str2num(exception{I}(length(exception{I})-3:length(exception{I})))];
    end
    temp=[];
    for I=1:size(y,1)
        if(isempty(find(numexc==I)))
            temp=[temp;y(I)];
        end
    end
    y=temp;
end
%-----------------------
for k=2:5
    warning off;
    DirName=strcat(Keys.root_package,'\clustering\',path);
    mkdir(DirName);
    warning on;
    DirName=strcat(DirName,'\',num2str(size(dir(strcat(Keys.root_package,'\clustering\',path)),1)-1));
    mkdir(DirName);
    IDX = kmeans(matrix,k);
    for J=1:k
        cluster=matrix(find(IDX==J),:);
        cluster_act=y(find(IDX==J));
        save(strcat(DirName,'\cluster',num2str(J),'.mat'),'cluster');
        save(strcat(DirName,'\cluster_act',num2str(J),'.mat'),'cluster_act');
    end
    method_param=strrep(strcat('kmeans("',Keys.root_package,'\MDmatrix\',path,'\matrix.mat",', num2str(k),')'),'\','\\');
    fid=fopen(strcat(DirName,'\method_param.txt'),'w');
    fprintf(fid,method_param);
    fclose(fid);
   
end

end