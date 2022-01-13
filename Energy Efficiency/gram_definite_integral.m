function Wm = gram_definite_integral(A,B,T)
%   Calculate gram matrix by definite integral
%   Slow but also works on uncontrollable system
%
%   W = \int_{0}^{T} e^{A t}BB^Te^{A^T t}dt
%
%   Inputs:     A,     N-by-N matrix, system matrix
%               B,     N-by-M matrix, control matrix with M control nodes
%               T,     Constant, control time
%
%   Output:     Wm,    N-by-M matrix, system gram matrix
    fun = @(x) expm(A*x)*(B*B')*expm(A'*x);
    Wm = integral(fun,0,T,'ArrayValued',true);

end

