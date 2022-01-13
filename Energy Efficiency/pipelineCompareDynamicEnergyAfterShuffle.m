function Er = pipelineCompareDynamicEnergyAfterShuffle(time_series, control_node, window_size, r_time)
%   Pipeline of Compare dynamic control energy before and after shuffle
%    dynamicMiniEnergy to calculate the control energy.
%   This pipeline include:
%       1. random select the start point and final point
%       2. calculate the snapshot number and dynamic functional connectivities
%       3. normalize the functional connectivities to system matrix
%       4. calculate the dynamic minimum control energy.
%       5. random shuffle the snapshot order and recompute the mimimum control
%       energy
%   Inputs:     time_series,     N-by-T-by-M  matrix, fMRI BOLD series, N
%                               regions, T time points and M subjects
%               control_node,     N-by-1 vector, control number, e.g. [2, 5, 7, 10]
%               delta_t,     Constant, control time for one snapshot
%               x0,     N-by-1 vector, which means the initial state
%               xf,     N-by-1 vector, which means the final state
%               r_time  Constant, repeat times
%
%   Output:     Er,      (r_time+1)-by-1 vector, first one is minimum dynamic control
%                   energy and the rest is minimum dynamic control energy
%                   after disrupt the order of sanpshots
%   Reference: Li, A., Cornelius, S. P., et, al.  Science (2017) 358, 1042¨C1046
    [N,T,M] = size(time_series);
    B = eye(N);
    B = B(:,control_node);
    Er = zeros(M,r_time+1);
    start_point = randi(floor(T/2.5)); % random select a start point
    num_snapshot = randi(floor((T-start_point-window_size)/window_size)); % random choose the number of snapshot
    final_point = start_point + window_size * num_snapshot;
    delta_s = randi(6000)/10; % random select a control time for static control
    delta_t = delta_s / num_snapshot; % dynamic control time for each snapshot
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
        
        % calculate the dynamic functional connectivity
        
        dFCM = zeros(N,N,num_snapshot);
        for di = 1:num_snapshot
            [dFCi,p] = corr(Si(:,start_point+window_size*(di-1):start_point+window_size*di)');
            dFCi(p>0.05) = 0;
            dFCM(:,:,di) = dFCi;
        end
        
        % system matrix
        dAM = functional2system(dFCM);
        
        Er(sbn,:) = dynamicMiniEnergyDisruptOrder(dAM,B,delta_t,x0,xf,r_time);
    end
end

