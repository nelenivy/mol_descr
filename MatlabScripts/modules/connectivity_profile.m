function connectivity_profile(KeyFile)
% sozdaet profile dlya viborki
% -ALL.prf - vse moleculy viborki
% SingleConnectivity.prf - vse odnosvyaznie moleculy

Keys=readKeys(KeyFile);
exception={};
warning on;
curDir=pwd;
if (isempty(Keys.profile_file))
    profile=strcat(curDir,'\',Keys.root_package,'\molecules\',Keys.root_package, '-', Keys.mol_prefix, '-ALL.prf' )
    PROFILEid=fopen(profile,'w');
else
    profile=strcat(curDir,'\',Keys.root_package,'\molecules\', Keys.root_package,'-', Keys.mol_prefix, '-', Keys.profile_file, '.prf' );
    NotSingleConnectivity=strcat(curDir,'\',Keys.root_package,'\molecules\', Keys.root_package,'-', Keys.mol_prefix, '-', 'NotSingleConnectivity.prf' );
    PROFILEid=fopen(profile, 'w');
    NotSingleConnectivityID=fopen(NotSingleConnectivity, 'w');
end


fprintf(PROFILEid, strcat('root_package=', Keys.root_package,'\n', 'mol_prefix=',Keys.mol_prefix,'\n','profile_file=',Keys.profile_file,'\n'));

fprintf(NotSingleConnectivityID, strcat('root_package=',Keys.root_package,'\n','mol_prefix=',Keys.mol_prefix,'\n','profile_file=NotSingleConnectivity \n'));

files=dir(strcat(curDir,'\', Keys.root_package,'\molecules'));
n=size(files);
for I=3:n(1)
    
    k=findstr(files(I).name, Keys.mol_prefix);
    if (~isempty(k))
        if ((k(1)==1) && (files(I).isdir==1))
            fileName=strcat(curDir,'\', Keys.root_package,'\molecules\', files(I).name,'\',files(I).name, '.mol');
            fid_n=fopen(fileName,'r');
            if (fid_n)
                CMatrix=ConnectMatrix(fileName);
                components=ConComponents(CMatrix);
                if (components==1)
                    fprintf(PROFILEid, files(I).name);
                    fprintf(PROFILEid, '\n');
                else
                    warning (strcat('molecula  "', files(I).name, '" is not 1-connected.'));
%                    choice=input('Press the key: \n 1. Exclude the molecule from the consideration \n 2. Continue, ignore the nonsimple connectivity \n');
                    choice=2;
                    if (choice==1)
                        exception=cat(2,exception, files(I).name);
                    end
                    fprintf(NotSingleConnectivityID, files(I).name);
                    fprintf(NotSingleConnectivityID, '\n');
                end
            end
            fclose(fid_n);
        end
    end
end

save(strcat(curDir,'\', Keys.root_package,'\exceptions.mat'),'exception');

fclose('all');
