function Wm = gram_lyaplov(A,B,T)
%   Calculate gram matrix by Lyaplov equation (controllable system)
%
%   Solve Lyaplov equation: AW + WA^T +BB^T = 0
%
%   Inputs:     A,     N-by-N matrix, system matrix
%               B,     N-by-M matrix, control matrix with M control nodes
%               T,     Constant, control time
%
%   Output:     Wm,    N-by-M matrix, system gram matrix
    [row,col] = size(A);
    C = eye(row,col);
    D = zeros(size(B));
    sys =  ss(A,B,C,D); % Build a linear systems
    opt = gramOptions('TimeInterVals',[0,T]); % Control during finite time
    Wm = gram(sys,'c',opt);
end

