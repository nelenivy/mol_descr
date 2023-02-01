function Descriptor_all(KeyFile)
%vyzyvaet metod Descriptor dlya vseh faylov, raspologennih v kornevoy papke
% root_package i prefiksom mol_prefix. Vhod - fayl s klyuchami
Keys=readKeys(KeyFile);

% proverit Keys
if((isempty(Keys.root_package)) || (isempty(Keys.mol_prefix)) ||(isempty(Keys.chain_length))|| (isempty(Keys.markers)) || (isempty(Keys.distance_type)))
    error('The Key-file contains not enough data');
end
load(strcat(pwd,'\',Keys.root_package,'\exceptions.mat'));

files=dir(strcat(pwd,'\', Keys.root_package,'\molecules')); %imena fseh direktoriy
n=size(files);
for I=3:n(1)
    k=findstr(files(I).name, Keys.mol_prefix);
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1))
            NewDir=strcat(pwd,'\',Keys.root_package,'\molecules\',files(I).name);
            fileName=strcat(NewDir,'\',files(I).name, '.mol'); % formiruem im9 nujnogo fayla
            %fileName
            if (fopen(fileName,'r') && isempty(find(strcmp(exception,files(I).name))))
               Descriptor(fileName, str2num(Keys.chain_length), Keys.markers);
            end
            fclose('all');
        end
    end
end

clear;
fclose('all');