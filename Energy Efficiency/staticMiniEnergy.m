function E = staticMiniEnergy(As,B,T,x0,xf)
%   Calculate the minimum static control energy under unrestricted control path
% 
%   E(x_0, x_f) = \frac{1}{2} d^T W^{-1}d 
% 
%   Inputs:     As,     N-by-N matrix, system matrix
%               B,     N-by-M matrix, control matrix with M control nodes
%               T,     Constant, control time
%               x0,     N-by-1 vector, which means the initial state
%               xf,     N-by-1 vector, which means the final state
%
%   Output:     E,      Constant, minimum static control energy
%   Reference: Li, A., Cornelius, S. P., et, al.  Science (2017) 358, 1042¨C1046
    try
        Ws = gram_lyaplov(As,B,T); % calculate gram matrix by Lyaplov equation (controllable system)
    catch
        Ws = gram_definite_integral(As,B,T); % calculate gram matrix by definite integral
    end
    d = xf - expm(As*T)*x0; % distance between the final state and the desired state
    E = 0.5*d'*pinv(Ws)*d;

end

