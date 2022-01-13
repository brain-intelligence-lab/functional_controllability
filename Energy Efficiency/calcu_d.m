function d = calcu_d(EA, x0, xf)
%   Calculate the distance between the final point 
%   and the desired point without input energy.
% 
%   d = x_f - \prod_{i=M}^{1} e^(A_i \Delta t_i) x_0
% 
%   Inputs:     EA,     N-by-N-by-M matrix, where contains M matrix e^(A_i \Delta t_i)
%               x0,     N-by-1 vector, which means the start state
%               xf,     N-by-1 vector, which means the final state
%
%   Output:     d,      N-by-1 vector, state distance
    [row,col,num] = size(EA);
    M = eye(row,col);
    for i = 1:num
        M = EA(:,:,i)*M; % \prod_{i=M}^{1} e^(A_i \Delta t_i)
    end
    d = xf-M*x0;

end

