function MDmatrix(KeyFile, descriptor_list)
%создает MD матрицу
display('MDmatrix...');
Keys=readKeys(KeyFile);
% проверить Keys
if ((isempty(Keys.root_package)) || (isempty(Keys.mol_prefix)) || (isempty(Keys.chain_length) || isempty(Keys.method)))
    return;
end
MD=[];

    load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));
   
if (nargin==1)
    load(strcat(pwd,'\',Keys.root_package,'\descriptor_list\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,'\descriptor_list.mat'));
end
    n=size(descriptor_list);

files=dir(strcat(Keys.root_package,'\molecules'));

for I=3:size(files,1)
    k=findstr(files(I).name, Keys.mol_prefix);
    
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1) && isempty(find(strcmp(exception,files(I).name))))
                  
            % имя открываемого файла
            load(strcat(pwd,'\',Keys.root_package,'\molecules\',files(I).name,'\',Keys.method,'\',Keys.chain_length,'\', Keys.markers,'\descriptors.mat'));
            MD=[MD; zeros(1, n(1))];
            
            %заполнение I строки матрицы
            for K=1:size(descriptors,2)
                ind=binary_search(descriptor_list, descriptors(K).chain);% ~=0
                if (ind >0)
                    MD(size(MD,1), ind)=descriptors(K).count; 
                end
            end
            clear descriptors
        end
    end
end
matrix=MD;
warning off;
mkdir(strcat(Keys.root_package,'\MDmatrix\',Keys.method,'\', Keys.chain_length,'\',Keys.markers));
warning on;
save(strcat(Keys.root_package,'\MDmatrix\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,'\matrix.mat'),'matrix');
fclose('all');
