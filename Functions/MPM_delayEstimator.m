function vector_num_estimated_paths = MPM_delayEstimator(multipath_channel_all_packets,num_packets,antennas,P,unit_circle_toll)
vector_num_estimated_paths = NaN(1,num_packets);
for num_packet = 1:num_packets
        
    all_gains = NaN.*ones(antennas, 100);
    num_estimated_paths = NaN(1,antennas);
        
    H_iscell = iscell(multipath_channel_all_packets);
    if (H_iscell)
        H_tmp = multipath_channel_all_packets{num_packet};
    else
        H_tmp = squeeze(multipath_channel_all_packets(num_packet, :, :, :));
    end
    
    multipath_channel = H_tmp(:, :, 1);
    num_subcarriers = size(multipath_channel,1);

    for antenna = 1:antennas
        H_n = multipath_channel(:,antenna);
        
        % Hankel matrix:
        K = num_subcarriers-P;
        Hankel_matrix = [];
        for k = 1:K
            Hankel_matrix(k,:) = H_n(k:k+P);
        end
        H1 = Hankel_matrix(:,1:P);
        H2 = Hankel_matrix(:,2:P+1);
        H_H = pinv(H2)*H1;
        z_l = eig(H_H);
        
        estimated_delay = (num_subcarriers/(2*pi)).*(angle(z_l));
        
        %----distance to the unit circle----%
        d = abs(1-abs(z_l));
        %-----------------------------------%
        
        ind_d = find(d<unit_circle_toll);
        if isempty(ind_d)
            [~, ind_d] = min(d);
        end
        delays = sort(estimated_delay(ind_d));
        A_tau = [];
        for n = 1:num_subcarriers
            z_l_per_subcarrier = zeros(1, length(delays));
            for l = 1:length(delays)
                z_l_per_subcarrier(l) = exp(-1i*2*pi*delays(l)*n/num_subcarriers);
            end
            A_tau = vertcat(A_tau, z_l_per_subcarrier);
        end        
        gains = pinv(A_tau) * H_n;
        all_gains(antenna,1:length(gains)) = transpose(gains);
        num_estimated_paths(antenna) = length(gains);
    end
    vector_num_estimated_paths(num_packet) = floor(mode(num_estimated_paths));
end


