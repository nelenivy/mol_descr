function cluster_clusterdata(KeyFile, YFile)

% функци€ кластеризации матрицы методом иерархического кластерного анализа. –езультат сохран€етс€ в
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

% обработка несв€зных молекул
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

distance={'euclidean','seuclidean','cityblock','minkowski','hamming','jaccard','chebychev'};
linkage={'average','centroid','complete','median','single','ward','weighted'};
method_param=[];
for I1=1:length(distance)
    for I2=1:length(linkage)
        
        DirName=strcat(Keys.root_package,'\clustering\',path);
        mkdir(DirName);
        
        DirName=strcat(DirName,'\',num2str(size(dir(strcat(Keys.root_package,'\clustering\',path)),1)-1));
        mkdir(DirName);
        
        T=clusterdata(matrix,'distance',distance{I1},'linkage',linkage{I2},'maxclust',5);
        
        for J=1:5
            cluster=matrix(find(T==J),:);
            cluster_act=y(find(T==J));
            save(strcat(DirName,'\cluster',num2str(J),'.mat'),'cluster');
            save(strcat(DirName,'\cluster_act',num2str(J),'.mat'),'cluster_act');
        end
        method_param=strrep(strcat('clusterdata(',Keys.root_package,'\MDmatrix\',path,'\matrix.mat,','distance,',distance{I1},',linkage,',linkage{I2},',maxclust,',num2str(5),')'),'\','\\');
        fid=fopen(strcat(DirName,'\method_param.txt'),'w');
        fprintf(fid,method_param);
        fclose(fid);
        %save(strcat(DirName,'\method_param.mat'),'method_param');
    end
end