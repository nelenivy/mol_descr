function res=distance(Mol, cluster)

% ������� ������� ���������� �� �������� �� ������ ��������
% ������� ���������
if (size(cluster,1)>1)
    cluster=mean(cluster);
end
res=0;
for I=1:size(Mol,2)
    res=res+(Mol(I)-cluster(I))^2;
end
res=sqrt(res);
