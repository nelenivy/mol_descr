function OT2main(KeyFile, DistFile, SDFFile)
%Distance=
%метод дл€ построени€ особых точек второго уровн€, заупскаетс€ после того,
%как построены ќ“ 1 уровн€

%вход - файл с параметрами запуска KeyFile
%DistFile - значени€ и имена интервалов, рассто€ние между цепочками
display('OT2main...');

if nargin<3
    SDFFile='0';
end
%проверить наличие входных параметров
if (~exist(KeyFile))
    error('Please, specify the path to KeyFile');
end
if (~exist(DistFile))
    error('Please,specify the path to DistanceFile');
end

% если не выборка не обрабатывалась, запустить обработку с нул€
Keys=readKeys(KeyFile);
if (~strcmp(SDFFile,'0') && ~(isdir(strcat(Keys.root_package,'\molecules'))))
    cutSDF(SDFFile, KeyFile);
else
    if (~isdir(strcat(Keys.root_package)))
        error('Training sample should be submitted as a sdf-file or directory in the root');
    end
end
if (~exist(strcat(Keys.root_package,'\exceptions.mat')))
    connectivity_profile(KeyFile);
end


Distance=OT2_all(KeyFile, DistFile);
out=CodeAllOT(KeyFile, DistFile);
if (out==0)
    Descriptor_list(KeyFile);
    MDmatrix(KeyFile);
end

end