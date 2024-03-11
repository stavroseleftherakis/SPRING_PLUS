function [median_MUSIC_AoAestimation,FirstPath_AoAestimation,MUSIC_2DSpectrum_AllPcks] = get_AoAestimations(H,num_packets,deltaPhi,SubCarrInd,T,c,fc,fgap,space_between_antennas,SpotFi_structure,NumEstPaths_MPM)

vector_MUSIC_AoAestimations = NaN(1,num_packets);

for num_packet = 1:num_packets

    %----CSI per packet----%
    H_iscell = iscell(H);
    if (H_iscell)
        H_tmp = H{num_packet};
    else
        H_tmp = squeeze(H(num_packet, :, :, :));
    end
    scaled_H = H_tmp(:, :, 1);
    CSI = transpose(scaled_H);
    %-----------------------%
    
    %----SPRING calibration----%
    calibratedCSI = get_calibratedCSI(CSI,deltaPhi);
    %--------------------------%
    
    % ----Algorithms for AoA estimations----
    %|                                      |
    %|       Normal MUSIC estimator         |
    R = cov(calibratedCSI');
    [AoA_Estimation,spec,specang] = musicdoa(R,1,'ScanAngles',-90:0.1:90);
    vector_MUSIC_AoAestimations(num_packet) = AoA_Estimation;
    %|                                      |
    %|     1st-FirstPath AoA Estimator      |
    [matrix_3Dinfo_SpotFi_perPck(:,:,num_packet), MUSIC_2DSpectrum_AllPcks(:,:,num_packet)] = get_matrix_3Dinfo_SpotFi_PathAdaptation(calibratedCSI,SubCarrInd,T,c,fc,fgap,space_between_antennas,SpotFi_structure,NumEstPaths_MPM(num_packet));
    %|                                      |
    % --------------------------------------
    
end

%----Storing MUSIC----%
vector_MUSIC_AoAestimations(isnan(vector_MUSIC_AoAestimations)) = [];
median_MUSIC_AoAestimation = median(vector_MUSIC_AoAestimations);
%----------------------------%

%----2nd-FirstPath AoA Estimator----%
matrix_3Dinfo_SpotFi(:,:) = matrix_3Dinfo_SpotFi_perPck(:,:,:);
%-----------------------------------%

%-----------------------------------%
FirstPath_AoAestimation = FirstPath_AoAestimator(matrix_3Dinfo_SpotFi,SpotFi_structure);
%-----------------------------------%