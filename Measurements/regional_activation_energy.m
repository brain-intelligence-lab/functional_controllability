function E = regional_activation_energy(FC)
% Compute the regional activation energy for every regions under whole
% brain control.
% This function need to import the function in
% functional_controllability/Energy Efficiency
%
%   inputs:     FC,     N-by-N, connectivity matrix
%
%   outputs:    E,      N-by-1 vector, static minimum control energy for
%   corresponding regions activation x0=0 to xf=1
    A = functional2system(FC);
    N = size(A,1);
    E = zeros(N,1);
    B = eye(N);
    x0 = zeros(N,1);
    T = 1;
    for i = 1: N
        xf = zeros(N,1);
        xf(i) = 1;
        E(i) = staticMiniEnergy(A,B,T,x0,xf);
    end
end

