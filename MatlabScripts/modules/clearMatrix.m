function [newMatrix, newY, newR2, newcoef, newInf,newexcluded]=clearMatrix(Matrix, Y, R2, coef, inf,excluded)

%������� �� ��������� ������ ������

%����:
%Matrxi - ��-�������
%R2 - �������� ��������
%coef - ������������ �������������� �������
%inf - 

%�����:
%������� ��������� ��� ������ �����

newMatrix={};
newY={};
newR2={};
newcoef={};
newInf={};
newexcluded={};
for I=1:size(Matrix,2)
    if (~isempty(Matrix{I}))
       newMatrix=cat(2,newMatrix,Matrix{I});
       newY=cat(2, newY,Y{I});
       newR2=cat(2, newR2, R2{I});
       newcoef=cat(2, newcoef, coef{I});
       newInf=cat(2, newInf, inf{I});
       newexcluded=cat(2,newexcluded,excluded{I});
    end
end