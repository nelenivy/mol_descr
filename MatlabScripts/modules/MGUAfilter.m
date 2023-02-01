function [reducedMatrix,excluded]=MGUAfilter(matrix)
% функция предобработки матрицы перед запуском МГУА
% если в столбце количество ненулевых значений меньше 5%, складываем
% значение этого столбца с первым столбцом результирующей матрицы
% номера неинформативных столбцов в excluded

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
