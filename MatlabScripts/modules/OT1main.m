function OT1main(KeyFile, SDFFile)
%
%  ���������� ��-�������, ����������� ������� ������, ��� ������� ��������
%  � ��������� SDFFile
%  ��������� ��������� �������� � ��������� ����� KeyFile
%  ������ KeyFile:
%       root_package=BETUL
%       mol_prefix=betul
%       profile-file=SingleConnectivity
%       Chain_length=3
%       distance_type=geometric
%       markers=dbr
% 
%root_package - ��� ����������, � ������� ��������� ���������� mol-�����,
% � ����� ������ ������������ � ��-�������
%mol_prefix - ������ ����� ������������ mol-�����
%profile-file - �������� �������. � ������ ������ - ��������� ���������
%Chain_length - ����� ����������� ������� ������� ������. ��������: 2, 3, 4
%distance_type - ��� ������������ ����������. ��������: geometric/topologic
%markers - ��������� ��������. ��������: ___, __r, _b_, d__, db_, d_r, _br, dbr


% ����:
% KeyFile - ��������� ���� � ����������� ���������
% SDFFile - �������������� �������� - ��� ����� sdf � ��������������
% �������� (���� ������� ��� ������������ � ���� ������� ������ ����������, �� ��������� �� ����������, ��������� � KeyFile)
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
   