function Wd = diagW(WM)
% build huge diag gram matrix
%
% W = diag(W_1, W_2, ..., W_M)
%
%   Inputs:     WM,     N-by-N-by-M matrix, M N-by-N gram matrix
%
%   Output:     Wd,      N*M-by-N*M diag gram matrix
    [row,col,num] = size(WM);
    B = zeros(row,col);
    P = cell(num);
    for i = 1:num
        for j = 1:num
            if i==j
                P(i,j)={WM(:,:,i)};
            else
                P(i,j)={B};
            end
        end
    end
    Wd = cell2mat(P); 


end

