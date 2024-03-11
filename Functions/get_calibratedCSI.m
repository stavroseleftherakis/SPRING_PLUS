function calibratedCSI = get_calibratedCSI(CSI,delta_phi)

% Need description...

[num_antennas, num_subcarriers] = size(CSI); 
calibratedCSI = zeros(num_antennas, num_subcarriers);
total_mean = zeros(num_antennas,num_subcarriers);
single_total_mean = zeros(1,num_antennas);
for idx_antenna = 1:num_antennas
    CSI_per_antenna = CSI(idx_antenna,:); 
    CSI_frequency = fft(abs(CSI_per_antenna));
    W_f(1:8) = CSI_frequency(1:8);
    W_f(12:224) = CSI_frequency(12:224);
    W_f(228:234) = CSI_frequency(228:234);
    calibrated_amplitude_CSI = abs(ifft(W_f));
    single_total_mean(idx_antenna) = mean(calibrated_amplitude_CSI);
    total_mean(idx_antenna,:) = ones(1,num_subcarriers)*mean(calibrated_amplitude_CSI);
    calibrated_CSI_per_antenna = calibrated_amplitude_CSI.*exp(1i*angle(CSI_per_antenna)).*exp(-1i*delta_phi(idx_antenna));
    calibratedCSI(idx_antenna,:) = calibrated_CSI_per_antenna;
end

calibratedCSI = mean(single_total_mean)*calibratedCSI./total_mean;

end
