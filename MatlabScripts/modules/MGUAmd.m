function MGUAmd(KeyFile, YFile)

% метод МГУА для МД-матриц, параметры и путь к файлу хранятся в KeyFile
display('MGUAmd...');
Keys=readKeys(KeyFile);
load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));
warning off;
path=strcat(Keys.method,'\',Keys.chain_length,'\',Keys.markers)
fid=fopen(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'),'r');
if(fid<0)
    error('MDMatrix is not exist');
end
fclose(fid);
load(strcat(Keys.root_package,'\MDmatrix\',path,'\matrix.mat'));
y=load(YFile);
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
y=y';
descNum=max(2,ceil(5*size(matrix,1)/100));
fid=fopen(strcat(Keys.root_package,'\MGUA\MDmatrix\',path,'\MGUA.mat'),'r');
if (fid<0)
    [reducedMatrix,excluded]=MGUAfilter(matrix);
    
    if (size(reducedMatrix,2)>descNum)
        descNum=max(2,ceil(5*size(reducedMatrix,1)/100));
        [result, used, coefficients, R2, informative] = f_mgua_v2(reducedMatrix, y, descNum, descNum, 0.005, 0.99, 0.85, 1);
    else
        excluded=0;
        [result, used, coefficients, R2, informative] = f_mgua_v2(matrix, y, descNum, descNum, 0.005, 0.99, 0.85, 1);
    end
warning off;
mkdir(strcat(Keys.root_package,'\MGUA\MDmatrix\',path));
warning on;
save(strcat(Keys.root_package,'\MGUA\MDmatrix\',path,'\MGUA.mat'),'coefficients','R2','informative','excluded');
else
    fclose(fid);
end

end