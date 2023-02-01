function generalCluster(KeyFile,param)

thresholdR2=0.01; % порог для R2, если R2 кластеров отличаются менее, чем на это значение, то не размельчаем
thresholdSize=3; % порог для размера кластера, если полученный кластер меньше этого размера, то не выделяем в отдельный кластер
display('Clustering...');
Keys=readKeys(KeyFile);
clusterMatr={}; % матрицы кластеров
clusterR2={}; % R2 соответствующего кластера
clusterCoef={}; % коэффициенты лучшей модели для каждого кластера
clusterInf={};
clusterY={}; % вектор значений активности для кластера
clusterMethodInd={}; % индекс метода, в котором получен лучший результат
clusterExcluded={};
    %для кластера
path=strcat('\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,'\');
    
    % инициализация кластерных данных значениями первого метода
    % кластеризации
N=0;
while(exist(strcat(Keys.root_package,'\clustering',path ,'1\cluster',num2str(N+1),'.mat'))==2)
     N=N+1;
end
    
for J=1:N
    cluster1=load(strcat(Keys.root_package,'\clustering',path ,'1\cluster',num2str(J),'.mat')); %матрица кластера
    load(strcat(Keys.root_package,'\clustering',path ,'1\cluster_act',num2str(J),'.mat')); % вектор активности для кластера
    clusterMatr{J}=cluster1.cluster;
    clusterY{J}=cluster_act;
    load(strcat(Keys.root_package,'\MGUA\clustering',path,'1\MGUA',num2str(J),'.mat')); % результаты МГУА для кластера
    clusterR2{J}=R2;
    clusterCoef{J}=coefficients;
    clusterInf{J}=informative;
    clusterExcluded{J}=excluded;
    clusterMethodInd{J}='1';
        
end

for J=2:size(dir(strcat(Keys.root_package,'\clustering',path)),1)-2 % всего J методов кластеризации
    
    N=0; % количество кластеров при текущей кластеризации
    while(exist(strcat(Keys.root_package,'\clustering',path ,num2str(J),'\cluster',num2str(N+1),'.mat'))==2)
        N=N+1;
    end
        % рассматриваем все возможные пары кластеров, ищем пересечения
    for J1=1:N
        cluster1=load(strcat(Keys.root_package,'\clustering',path ,num2str(J),'\cluster',num2str(J1),'.mat'));
        curCluster=cluster1.cluster;
        load(strcat(Keys.root_package,'\clustering',path ,num2str(J),'\cluster_act',num2str(J1),'.mat'));
        curY=cluster_act;
        load(strcat(Keys.root_package,'\MGUA\clustering',path,num2str(J),'\MGUA',num2str(J1),'.mat'));
        reserve=clusterMatr;
        reserveY=clusterY;
        for J2=1:size(clusterMatr,2)
                % сравнить R2 нового кластера с R2 текущих кластеров
                % пересечение кластеров принадлежит кластеру с большим R2
            if (isempty(R2))
                R2=0;
            end
            if (R2>(clusterR2{J2}+thresholdR2) && ~isempty(curCluster))
                [curCluster,clusterMatr{J2},clusterY{J2}]=intersection(curCluster, clusterMatr{J2},clusterY{J2});
            else
                [clusterMatr{J2},curCluster,curY]=intersection(clusterMatr{J2},curCluster,curY);
            end
                
         end
            if (size(clusterMatr,2)==5)
                a='sdf';
            end
         if (~isempty(curCluster) && size(curCluster,1)>thresholdSize) % сохранить новый кластер
             [clusterMatr,clusterY, clusterR2, clusterCoef,clusterInf, clusterExcluded]=clearMatrix(clusterMatr,clusterY, clusterR2, clusterCoef, clusterInf,clusterExcluded);
             clusterMatr=cat(2,clusterMatr,curCluster);
             clusterY=cat(2,clusterY,curY);
             clusterR2=cat(2,clusterR2,R2);
             clusterCoef=cat(2,clusterCoef,coefficients);
             clusterInf=cat(2, clusterInf, informative);
             clusterExcluded=cat(2,clusterExcluded, excluded);
             clusterMethodInd=cat(2,clusterMethodInd, num2str(J));
          else
             clusterMatr=reserve;
             clusterY=reserveY;
          end
      end
end
 
 % сохранить результат
 mkdir(strcat(Keys.root_package,'\clustering',path,'generalClustering'));
 mkdir(strcat(Keys.root_package,'\MGUA\clustering',path,'generalClustering'));
 for J=1:size(clusterMatr,2)
     cluster=clusterMatr{J};
     save(strcat(Keys.root_package,'\clustering',path,'generalClustering\cluster',num2str(J),'.mat'),'cluster');
     cluster_act=clusterY{J};
     save(strcat(Keys.root_package,'\clustering',path,'generalClustering\cluster_act',num2str(J),'.mat'),'cluster_act');
     R2=clusterR2{J};
     informative=clusterInf{J};
     coefficients=clusterCoef{J};
     excluded=clusterExcluded{J};
     save(strcat(Keys.root_package,'\MGUA\clustering',path,'generalClustering\MGUA',num2str(J),'.mat'),'coefficients','R2','informative','excluded')
 end