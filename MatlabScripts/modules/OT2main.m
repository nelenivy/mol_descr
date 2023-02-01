function OT2main(KeyFile, DistFile, SDFFile)
%Distance=
%����� ��� ���������� ������ ����� ������� ������, ����������� ����� ����,
%��� ��������� �� 1 ������

%���� - ���� � ����������� ������� KeyFile
%DistFile - �������� � ����� ����������, ���������� ����� ���������
display('OT2main...');

if nargin<3
    SDFFile='0';
end
%��������� ������� ������� ����������
if (~exist(KeyFile))
    error('Please, specify the path to KeyFile');
end
if (~exist(DistFile))
    error('Please,specify the path to DistanceFile');
end

% ���� �� ������� �� ��������������, ��������� ��������� � ����
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