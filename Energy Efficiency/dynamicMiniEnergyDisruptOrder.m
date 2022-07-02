function E = dynamicMiniEnergyDisruptOrder(A,B,delta_t,x0,xf,r_time)
%   Calculate the minimum dynamic control energy before and after disrupt
%   the snapshort order
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
%               r_time  Constant, repeat times
%
%   Output:     E,      (r_time+1)-by-1 vector, first one is minimum dynamic control
%                   energy and the rest is minimum dynamic control energy
%                   after disrupt the order of sanpshots
    [row,col,num] = size(A);
    Wm = zeros(row,col,num);
    E = zeros(1,r_time+1);
    EA = zeros(row,col,num);  
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
    E(1) = 0.5*d'*pinv(Weff)*d; % dynamic minimum control energy befor shuffle
    
    % start shuffle r_times
    for ri = 1:r_time
        r_order = randperm(num);
        EAr = EA(:,:,r_order);
        Wmr = Wm(:,:,r_order);
        dr = calcu_d(EAr,x0,xf);
        Wr = diagW(Wmr);
        Sr = calcu_S(EAr);
        Weffr = Sr*Wr*Sr';
        E(ri+1) = 0.5*dr'*pinv(Weffr)*dr; % dynamic minimum control energy after shuffle
    end
end

