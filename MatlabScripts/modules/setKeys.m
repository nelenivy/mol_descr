function setKeys(keyFile, KeyStruct)
% � KeyFile ���������� �������� KeyStruct ��� ��������� ���������� �������

fid=fopen(keyFile,'w');
if (fid<0)
    error('Can not open file');
end
field=fields(KeyStruct);
for I=1:size(field,1)
    fprintf(fid, strrep(strcat(field{I},'=',KeyStruct.(field{I})),'\','\\'));
    fprintf(fid,'\n');
end

fclose(fid);