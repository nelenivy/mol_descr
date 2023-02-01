function Descriptor_list(KeyFile)
%vozvrawaet spisok vseh deskriptorov, vstrechaemih v vyborke
% zapisivaet v fayl mol_prefix.chK
% vhod - fayl s klyuchami
display('Descriptor_list...');
ChainMatrix={};
Keys=readKeys(KeyFile);
if (~isempty(Keys.root_package))
    curDir=strcat(pwd,'\',Keys.root_package,'\molecules');
   else return
end
if ((isempty(Keys.chain_length)) || (isempty(Keys.mol_prefix)) )
    return
end
    load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));
 
files=dir(curDir);
n=size(files);
for I=3:n(1)
    k=findstr(files(I).name, Keys.mol_prefix);
    %files(I).name
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1) && isempty(find(strcmp(exception,files(I).name))))
            NewDir=strcat(curDir, '\', files(I).name);
            fileName=strcat(NewDir,'\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,'\chains_coded.mat');  
            load(fileName);
            for J=1:size(chains_coded,1)
                % proverit est li takaya cepochka v matrice   
                res=binary_search(ChainMatrix, chains_coded(J,:));
                if (res==0)
                    ChainMatrix=cat(1, ChainMatrix, chains_coded(J,:));
                    ChainMatrix=sort(ChainMatrix);
                end                
            end 
        end
    end
   
end

descriptor_list=sort(ChainMatrix);
warning off;
mkdir(strcat(Keys.root_package,'\descriptor_list\',Keys.method,'\', Keys.chain_length,'\',Keys.markers));
warning on;
save(strcat(Keys.root_package,'\descriptor_list\',Keys.method,'\', Keys.chain_length,'\',Keys.markers,...
    '\descriptor_list.mat'),'descriptor_list');
  
fclose('all');