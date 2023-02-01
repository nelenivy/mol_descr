function total_activityCluster(KeyFile,param)

Keys=readKeys(KeyFile);
fid=fopen(strcat(Keys.root_package,'\predicted_activity\clustering\',Keys.method,'\total_activity.mat'),'r');
load(strcat(Keys.root_package,'\predicted_activity\clustering\',Keys.method,'\',Keys.chain_length,'\',Keys.markers,'\generalClustering\activity.mat'));
if (fid<0)
    total_activity=activity;
else
    load(strcat(Keys.root_package,'\predicted_activity\clustering\',Keys.method,'\total_activity.mat'));
    total_activity=total_activity+activity;
    fclose(fid);
end
save(strcat(Keys.root_package,'\predicted_activity\clustering\',Keys.method,'\total_activity.mat'),'total_activity');
