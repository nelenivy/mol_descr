function OT1main(KeyFile, SDFFile)
%
%  построение МД-матрицы, дескрипторы первого уровня, имя выборки задается
%  в параметре SDFFile
%  параметры обработки задаются в текстовом файле KeyFile
%  шаблон KeyFile:
%       root_package=BETUL
%       mol_prefix=betul
%       profile-file=SingleConnectivity
%       Chain_length=3
%       distance_type=geometric
%       markers=dbr
% 
%root_package - имя директории, в которой находятся полученные mol-файлы,
% а также список дескрипторов и МД-матрица
%mol_prefix - начало имени формируемого mol-файла
%profile-file - проверка условий. В данном случае - связность структуры
%Chain_length - длина формируемой цепочки связных атомов. Значения: 2, 3, 4
%distance_type - тип вычисляемого расстояния. Значения: geometric/topologic
%markers - включение маркеров. Значения: ___, __r, _b_, d__, db_, d_r, _br, dbr


% вход:
% KeyFile - текстовый файл с параметрами обработки
% SDFFile - необязательный параметр - имя файла sdf с обрабатываемой
% выборкой (если выборка уже представлена в виде нужного набора директорий, то обработка по параметрам, указанным в KeyFile)
%  
display('OT1.main...')
if nargin<2
    SDFFile='0';
end
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
Descriptor_all(KeyFile);
Descriptor_list(KeyFile);
MDmatrix(KeyFile);
clear;

end
   