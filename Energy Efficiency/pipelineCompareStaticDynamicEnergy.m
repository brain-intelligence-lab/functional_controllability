function [Es, Ed, delta_s] = pipelineCompareStaticDynamicEnergy(time_series, control_node, window_size)
%   Pipeline of Compare static and dynamic control energy
%   use staticMiniEnergy.m and dynamicMiniEnergy to calculate the control
%   energy.
%   This pipeline include:
%       1. random select the start point and final point
%       2. calculate the snapshot number and static and dynamic functional connectivities
%       3. normalize the functional connectivities to system matrix
%       4. calculate the static minimum control energy and dynamic minimum
%       control energy.
%   Inputs:     time_series,     N-by-T-by-M  matrix, fMRI BOLD series, N
%                               regions, T time points and M subjects
%               control_node,     N-by-1 vector, control number, e.g. [2, 5, 7, 10]
%               delta_t,     Constant, control time for one snapshot
%               x0,     N-by-1 vector, which means the initial state
%               xf,     N-by-1 vector, which means the final state
%
%   Output:     Es,      M-by-1 vector, minimum static control energy for
%                           subjects
%               Ed,      M-by-1 vector, minimum static control energy for
%                           subjects
%               delta_s,    control time \tau
%   Reference: Li, A., Cornelius, S. P., et, al.  Science (2017) 358, 1042¨C1046
    [N,T,M] = size(time_series);
    B = eye(N);
    B = B(:,control_node);
    Es = zeros(M,1);
    Ed = zeros(M,1);
    start_point = randi(floor(T/2.5)); % random select a start point
    num_snapshot = randi(floor((T-start_point-window_size)/window_size)); % random choose the number of snapshot
    final_point = start_point + window_size * num_snapshot;
    delta_s = randi(6000)/10; % random select a control time for static control
    delta_d = delta_s / num_snapshot; % dynamic control time for each snapshot
    
    for sbn = 1:M
        Si = time_series(:,:,sbn);
        Si = zscore(Si')'; % N-by-T
        % two ways to select the start point and final point
        % directly use the start point and final point
        x0 = Si(:,start_point);
        xf = Si(:,final_point);
        % calculate the mean value of the middle 10 point of the first and
        % last snapshots
%         x0 = mean(Si(:,start_point+window_size/2-10:start_point+window_size/2+10),2);
%         xf = mean(Si(:,final_point-window_size/2-10:final_point-window_size/2+10),2);
        
        % calculate the static and dynamic functional connectivity
        [sFC, p] = corr(Si(:,start_point:final_point)');
        sFC(p>0.05) = 0;
        
        dFCM = zeros(N,N,num_snapshot);
        for di = 1:num_snapshot
            [dFCi,p] = corr(Si(:,start_point+window_size*(di-1):start_point+window_size*di)');
            dFCi(p>0.05) = 0;
            dFCM(:,:,di) = dFCi;
        end
        
        % system matrix
        sA = functional2system(sFC);
        dAM = functional2system(dFCM);
        
        % control energy
        Es(sbn) = staticMiniEnergy(sA, B, delta_s, x0, xf);
        Ed(sbn) = dynamicMiniEnergy(dAM, B, delta_d, x0, xf);
        disp(['---finished subject:',num2str(sbn),' -----'])
end

