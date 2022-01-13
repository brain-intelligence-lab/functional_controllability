function SL = calcu_S(EA)
% S = (\prod_{i=M}^{2}e^{A_i\Dleta t},...,\prod_{i=M}^{t+1}e^{A_i\Dleta t},...I_N)
%   Inputs:     EA,     N-by-N matrix, e^{A_i\Dleta t}
%   Outputs:    SL,     N-by-N*M matrix
    [row,col,num] = size(EA);
    SL = eye(row,col);
    panel = eye(row,col);
    for k = 0:num-2
        p = num-k;
        panel = panel * (EA(:,:,p));
        SL = [panel SL];
    end
end

