function cutSDF(fileName, keys)
% po SDF faylu formiruet directoriyu, soderjawuyu otdelnie moleculy

line=' ';
counter=0;
OutDir='';
OutFile='';

fid=fopen(fileName, 'r');
Key=readKeys(keys);
if ((isempty(Key.root_package)) || (isempty(Key.mol_prefix)))
    return;
end
mkdir(Key.root_package);
mkdir(strcat(Key.root_package,'\molecules'));

curDir=pwd; % sozdanie profayla, soderjit vse moleculy
PROFILE='';
PROFILE=strcat(curDir,'\', Key.root_package,'\molecules\',Key.root_package, '-', Key.mol_prefix, '-ALL.prf')
PROFILEid=fopen(PROFILE, 'w');
fprintf(PROFILEid, strcat('root_package', '=',Key.root_package,'\nmol_prefix=',Key.mol_prefix,'\nprofile_file=ALL\n'));

while (line~=-1)
    counter=counter+1;
    % формирование директории
    OutDir=Key.mol_prefix;
    if ((counter>=1) && (counter<=9))
        OutDir=cat(2, OutDir, '000');
        OutDir=cat(2, OutDir, num2str(counter));
    end
    if ((counter>=10) && (counter<=99))
        OutDir=cat(2, OutDir, '00');
        OutDir=cat(2, OutDir, num2str(counter));
    end
    if ((counter>=100) && (counter<=999))
        OutDir=cat(2, OutDir, '0');
        OutDir=cat(2, OutDir, num2str(counter));
    end
    if ((counter>=1000) && (counter<=9999))
        OutDir=cat(2, OutDir, num2str(counter));
    end
    mkdir(strcat(curDir, '\',Key.root_package,'\molecules\', OutDir));
    
    %---------------------------------
    % создание нового файла
    OutFile='';
    OutFile=strcat(curDir,'\',Key.root_package,'\molecules\',OutDir, '\', OutDir, '.mol');
    OutId=fopen(OutFile, 'w');
        
    line=fgets(fid);
    %  mol-file
    while (line(1)~='$')
        if (line==-1)
            break
        end
        
        fprintf(OutId, line);
        line=fgets(fid);  
    end
    
    if (line~=-1)
        fprintf(PROFILEid, OutDir);
        fprintf(PROFILEid, '\n');
    end
    fclose(OutId);
end

fclose(PROFILEid);
rmdir(strcat(curDir,'\',Key.root_package,'\molecules\',OutDir),'s'); %udalit polednyuyu direktoriyu s pustym faylom
fclose('all');