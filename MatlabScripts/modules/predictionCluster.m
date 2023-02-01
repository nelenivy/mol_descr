function predictionCluster(KeyFile1, KeyFile2)

% предсказание активности для новых соединений
% путь к МД-матрицам новых соединений находится из KeyFile1
% результаты МГУА в KeyFile2

Keys1=readKeys(KeyFile1);
Keys2=readKeys(KeyFile2);
path=strcat('\',Keys1.method,'\',Keys1.chain_length,'\',Keys1.markers,'\');
N=0;
while(exist(strcat(Keys2.root_package,'\MGUA\clustering',path ,'generalClustering\MGUA',num2str(N+1),'.mat'))==2)
     N=N+1;
end
load(strcat(Keys1.root_package,'\MDmatrix',path,'matrix.mat'));
allR2={};
allCoefficients={};
allInformative={};
allExcluded={};
MolInClust=[];
activity=[];
R2Mol=[];
Dist2Clust=Inf(size(matrix,1), N);
for I=1:N
    cluster1=load(strcat(Keys2.root_package,'\clustering',path,'generalClustering\cluster',num2str(I),'.mat'));
    if (size(cluster1.cluster,1)>0)
        for M=1:size(matrix,1)
            Dist2Clust(M,I)=distance(matrix(M,:), cluster1.cluster);           
        end
        load(strcat(Keys2.root_package,'\MGUA\clustering',path,'generalClustering\MGUA',num2str(I),'.mat'));
        allR2=cat(1,allR2,R2);
        allCoefficients=cat(1,allCoefficients,coefficients);
        allInformative=cat(1,allInformative,informative);
        allExcluded=cat(1,allExcluded,excluded);
    end
end

for M=1:size(matrix,1)
    MolInClust=[MolInClust; find(Dist2Clust(M,:)==min(Dist2Clust(M,:)), 1)];
    ReducedMatrix=zeros(size(matrix,1),1);
    if (allExcluded{MolInClust(M)}~=0)
        for I=1:size(matrix,2)
            if (isempty(find(allExcluded{MolInClust(M)}==I)))
                ReducedMatrix=[ReducedMatrix, matrix(:,I)];
            else
                ReducedMatrix(:,1)=ReducedMatrix(:,1) + matrix(:,I);
            end
        end
    end
    R2Mol=[R2Mol; allR2(MolInClust(M))];
    activity=[activity; sign(allCoefficients{MolInClust(M)}* transpose([ReducedMatrix(M, allInformative{MolInClust(M)}), 1]))]; %вычисление активности
end
warning off;
mkdir(strcat(Keys1.root_package,'\predicted_activity\clustering',path,'generalClustering'));
warning on;
save(strcat(Keys1.root_package,'\predicted_activity\clustering',path,'generalClustering\activity.mat'),'activity');
