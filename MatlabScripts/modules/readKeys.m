function keyMatrix=readKeys(fileName)
keyMatrix=struct;
%('root_package','','mol_prefix','','distance_type','','Chain_length','','markers','','profile-file','');

fid=fopen(fileName, 'r');

line=fgets(fid);

while (line~=-1)
    name=strtrim(line(1:find(line=='=')-1));
    if(~isempty(name))
        line=strtrim(line(find(line=='=')+1:length(line)));
        if (~isempty(find(isspace(line),1)))
            keyMatrix.(name)=strtrim(line(1:find(isspace(line))));
        else
            keyMatrix.(name)=strtrim(line);
        end
        
    end
    line=fgets(fid);
end

if(~isfield(keyMatrix,'root_package'))
    error 'field "root_package" is empty';
elseif (~isfield(keyMatrix,'mol_prefix'))
    error 'field "mol_prefix" is empty';
elseif (~isfield(keyMatrix,'distance_type'))
    error 'field "distance_type" is empty';
elseif (~isfield(keyMatrix,'chain_length'))
    error 'field "Chain_length" is empty';
elseif(~isfield(keyMatrix,'markers'))
    error 'field "markers" is empty';
elseif (~isfield(keyMatrix,'profile_file'))
    error 'field "profile-file" is empty';
end
    

fclose(fid); 