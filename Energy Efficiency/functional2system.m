function A = functional2system(FC)
%   Translate functional connectivities to system matrix
%   A = -L/(\lambda_{max}(L))
%   L is the laplacian matrix for FC
%   Inputs:     FC,     N-by-N or N-by-N-by-M  matrix, functional matrix
%
%   Output:     A,      N-by-N, system matrix
    [row,~,num] = size(FC);
    if num == 1
        La = -FC;
        for i = 1:row
            La(i,i)= La(i,i) + sum(abs(La(i,:)));
        end
        A = -La / (max(abs(eig(La))));
    else
        La = -FC;
        A = zeros(row,row,num);
        for n = 1: num
            for i = 1:row
                La(i,i,n) = La(i,i,n) + sum(abs(La(i,:,n)));
            end
            A(:,:,n) = -La(:,:,n) / (max(abs(eig(La(:,:,n)))));
        end
    end
end

