function MGUAcluster(KeyFile,param)

% метод МГУА для матриц, полученных в результате кластеризации исходных
% матриц. Файлы хранятся в
% root_package\clustering\method\chain_length\markers\1(2,3,...)\clusterI.mat
% root_package\clustering\method\chain_length\markers\1(2,3,...)\cluster_actI.mat
display('MGUAcluster...');
Keys=readKeys(KeyFile);
warning off;
path=strcat(Keys.method,'\',Keys.chain_length,'\',Keys.markers);
for I=1:size(dir(strcat(Keys.root_package,'\clustering\',path)),1)
    for J=1:ceil(size(dir(strcat(Keys.root_package,'\clustering\',path,'\',num2str(I))),1)/2)
        fid=fopen(strcat(Keys.root_package,'\clustering\',path,'\',num2str(I),'\cluster',num2str(J),'.mat'),'r');
        if(fid>0)
            fclose(fid);
            warning off;
            %strcat(Keys.root_package,'\clustering\',path,'\',num2str(I),'\cluster',num2str(J),'.mat')
            cluster1=load(strcat(Keys.root_package,'\clustering\',path,'\',num2str(I),'\cluster',num2str(J),'.mat'));
            load(strcat(Keys.root_package,'\clustering\',path,'\',num2str(I),'\cluster_act',num2str(J),'.mat'));
            
                cluster_act=cluster_act';
                descNum=max(2,ceil(5*size(cluster1.cluster,1)/100));
                fid=fopen(strcat(Keys.root_package,'\MGUA\clustering\',path,'\',num2str(I),'\MGUA',num2str(J),'.mat'),'r');
                if (fid<0)
                    if (size(cluster1.cluster,1)>5 && length(find(cluster_act<0))>0 && length(find(cluster_act>0))>0)
                        [reducedMatrix,excluded]=MGUAfilter(cluster1.cluster);
                        
                        if (size(reducedMatrix,2)>=descNum)
                            descNum=max(2,ceil(5*size(reducedMatrix,1)/100));
                            [result, used, coefficients, R2, informative] = f_mgua_v2(reducedMatrix, cluster_act, descNum, descNum, 0.005, 0.99, 0.85, 1);
                        else
                            excluded=0;
                            [result, used, coefficients, R2, informative] = f_mgua_v2(cluster1.cluster, cluster_act, descNum, descNum, 0.005, 0.99, 0.85, 1);
                        end
                    else
                        result=[0]; used=[0]; coefficients=[0]; R2=[0]; informative=[];excluded=0;
                    end
                    warning off;
                    mkdir(strcat(Keys.root_package,'\MGUA\clustering\',path,'\',num2str(I)));
                    warning on;
                    save(strcat(Keys.root_package,'\MGUA\clustering\',path,'\',num2str(I),'\MGUA',num2str(J),'.mat'),'coefficients','R2','informative','excluded');
               else
                    fclose(fid);
               end
                          
            end
            
       
    end
end


end