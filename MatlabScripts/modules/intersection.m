function [cluster1, cluster3,Y3]=intersection(cluster1, cluster2,Y2)

% �������� 2 ��������, � ������� ������ �������� R2
% ������ - ����� ����������� � ������� ��� ����������� �� �������� �
% ������� R2
cluster3=[];
Y3=[];
for I=1:size(cluster2,1)
    if (IndOfMol(cluster2(I,:),cluster1)==0)
        cluster3=[cluster3; cluster2(I,:)];
        Y3=[Y3;Y2(I)];
    end
end