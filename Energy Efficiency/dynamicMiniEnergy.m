function E = dynamicMiniEnergy(A,B,delta_t,x0,xf)
%   Calculate the minimum dynamic control energy under unrestricted control path
%   Wm = \int_{0}^{T} e^{A_m t}BB^Te^{A_m^T t}dt
%   W = diag(W1, W2, ..., WM)
%   S = (\prod_{i=M}^{2}e^{A_i\Dleta t},...,
%       \prod_{i=M}^{t+1}e^{A_i\Dleta t},...I_N)
%   W_{eff} = SWS^T
%   d = x_f - \prod_{i=M}^{1} e^(A_i \Delta t_i) x_0
%   E(x_0, x_f) = \frac{1}{2} d^T W^{-1}d 
% 
%   Inputs:     A,     N-by-N-by-M  matrix, time-varying system matrix
%               B,     N-by-M matrix, control matrix with M control nodes
%               delta_t,     Constant, control time for one snapshot
%               x0,     N-by-1 vector, which means the initial state
%               xf,     N-by-1 vector, which means the final state
%
%   Output:     E,      Constant, minimum dynamic control energy
%   Reference: Li, A., Cornelius, S. P., et, al.  Science (2017) 358, 1042¨C1046

    [row,col,num] = size(A);
    Wm = zeros(row,col,num);
    EA = zeros(row,col,num); % e^{A_m t}
    for i = 1:num
        try
            Wm(:,:,i) = gram_lyaplov(A(:,:,i),B,delta_t);
        catch
            Wm(:,:,i) = gram_definite_integral(A(:,:,i),B,delta_t);
        end
        EA(:,:,i) = expm(A(:,:,i)*delta_t);
    end
    W = diagW(Wm);
    d = calcu_d(EA,x0,xf);
    S = longS(EA);
    Weff = S*W*S';
    E = 0.5*d'*pinv(Weff)*d;
end

