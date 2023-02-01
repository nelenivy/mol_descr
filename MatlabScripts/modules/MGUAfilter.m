function [reducedMatrix,excluded]=MGUAfilter(matrix)
% ������� ������������� ������� ����� �������� ����
% ���� � ������� ���������� ��������� �������� ������ 5%, ����������
% �������� ����� ������� � ������ �������� �������������� �������
% ������ ��������������� �������� � excluded

excluded=[];
reducedMatrix=zeros(size(matrix,1),1);
descNum=max(2,ceil(5*size(matrix,1)/100));
for I=1:size(matrix,2)
    NotNull=find(matrix(:,I)>0);
    if (length(NotNull)>descNum)
        reducedMatrix=[reducedMatrix, matrix(:,I)];
    else
        reducedMatrix(:,1)=reducedMatrix(:,1) + matrix(:,I);
        excluded=[excluded, I];
    end
end
