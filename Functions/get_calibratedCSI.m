function calibrated_CSI = get_calibratedCSI(CSI,delta_alpha_1,delta_alpha_2,delta_phi)
    [num_antennas, num_subcarriers] = size(CSI);
    calibrated_CSI = zeros(num_antennas, num_subcarriers);
    for idx_antenna = 1:num_antennas
        CSI_per_antenna = CSI(idx_antenna,:);
        calibrated_amplitude_CSI = db2mag((db(CSI_per_antenna) - delta_alpha_1(idx_antenna) - delta_alpha_2(idx_antenna,:)));
      
        calibrated_CSI_per_antenna = calibrated_amplitude_CSI.*exp(1i*angle(CSI_per_antenna)).*exp(-1i*delta_phi(idx_antenna));
        calibrated_CSI(idx_antenna,:) = calibrated_CSI_per_antenna;
    end
end
