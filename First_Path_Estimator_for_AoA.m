% --------------------------------------------------------------------%
%| This code implements SPRING+: 
%| S. Eleftherakis, G. Santaromita, M. Rea, X. Costa-Pérez, D. Giustiniano, 
%| "SPRING+: Smartphone Positioning from a Single WiFi Access Point", in 
%| IEEE Transactions on Mobile Computing, 2024.   
%| Corresponding author: Stavros Eleftherakis
%| Contact: stavros.eleftherakis@imdea.org
%|------------------------------------------------------------------- %

% cleaning environment
clear all
close all
clc
format
warning('off')

% Adding paths for important folders
addpath('Functions')
addpath('mat files')

% loading SPRING calibration mat files - This is needed only with SPRING+ HW (CSI_65_Examples.mat)
load('delta_alpha_1')
load('delta_alpha_2')
load('delta_phi')

% 802.11ac parameters
channels_802_11ac = [7:9 11:12 16:16:32 34:2:64 68 96 100:2:128 132:2:144 149:2:155];
frequency_channels_802_11ac = 5.035e9 +((channels_802_11ac-channels_802_11ac(1))*5)*10^6;
channel = 112;
fc = frequency_channels_802_11ac(channels_802_11ac==channel);
fgap = 312.5e3; % [Hz]
SubCarrInd = -122:122; % Indices of active subacarriers
index_null_subcarriers = [-103, -75, -39, -11, -1, 0, 1, 11, 39, 75, 103];
[~,idx] = intersect(SubCarrInd,index_null_subcarriers,'stable');
SubCarrInd(idx) = [];

% Fixed variables
c = .299792458e9; % speed of the light [m/s]
space_between_antennas = .026; % QTNA configuration

% Setting MPM inputs - The optimal values can vary depending on the scenario
MPM_PencilParameter = 5;
MPM_unitCircleToll = 0.35;

% Parameter to tune SF length
Nsl = 8;
K_tmp = [4 4 4 4 3 3 3 3]; % Hardcoded values. Can be automatized (Alg. 1)
J_tmp = [234 217 200 183 167 150 133 117]; % Hardcoded values. Can be automatized (Alg. 1)

input_data = load("mat files\CSI_65_Examples");

CSI_tmp = input_data.H; % CSI data
CSI = CSI_tmp(:,:,:,1); % SPRING+ HW have this behavior
[num_packets, num_subcarriers, num_antennas] = size(CSI);

SpotFi_structure = get_infoSpotFi; % SpotFi parameter inside the function
    
% MPM algorithm
estLpaths_for_PCK = MPM_delayEstimator(CSI,MPM_PencilParameter,MPM_unitCircleToll);
ToF_Difference_For_Consider_Same_Path = 10; % Sec. 4.3

% initialization 3D-matrix of estimation [path, [ToF, AoA, Corr], pck]
matrix_Estimated_perPck = NaN(SpotFi_structure.max_numPaths,3,num_packets);
    
for num_packet = 1:num_packets
    
    clear CSI_pkt
    
    CSI_pkt(:,:) = CSI(num_packet,:,:);
    % This calibration is needed only with SPRING+ HW (CSI_65_Examples.mat)
    CSI_pkt = get_calibratedCSI(CSI_pkt.',delta_alpha_1,delta_alpha_2,delta_phi);
    % CSI_pkt = CSI_pkt.' % If SPRING+ calibration is not used
    
    % Sanitization of SpotFi
    sanitized_and_in_line_CSI = SpotFi_Sanitization(CSI_pkt,num_subcarriers,num_antennas,SubCarrInd,1);
    
    maxNumPaths = SpotFi_structure.max_numPaths;
    
    % Alg. 1 implementation starts from HERE
    if num_packet == 1
        Index_SF_search = 1:Nsl; 
    else
        Index_SF_search = (Previus_r-1):(Previus_r+1);
        Index_SF_search((Index_SF_search < 1) | (Index_SF_search > 8)) = [];
    end
    
    for r = Index_SF_search
        SpotFi_structure.J = J_tmp(r);
        SpotFi_structure.K = K_tmp(r);
        [EstimatedMatrix,Spectrum] = backscatterEstimationMusic(sanitized_and_in_line_CSI.', num_antennas, num_subcarriers, c, fc,...
            SpotFi_structure.T, fgap, SubCarrInd, space_between_antennas, SpotFi_structure, ones(maxNumPaths));
        vector_corr = EstimatedMatrix(:,3);
        vector_corr = vector_corr./vector_corr(1); % corr normalization
        NumEstPath_SpotFi = sum(vector_corr >= SpotFi_structure.selected_corrThr);
        vector_ToF =  EstimatedMatrix(1:NumEstPath_SpotFi,1);
        if NumEstPath_SpotFi > 1
            if NumEstPath_SpotFi == 2
                if abs(vector_ToF(1) - vector_ToF(2)) <= ToF_Difference_For_Consider_Same_Path
                    NumEstPath_SpotFi = 1;
                end
            else
                counter = 0;
                for index_s = 1:(NumEstPath_SpotFi - 1)
                    if abs(vector_ToF(index_s) - vector_ToF(index_s+1)) <= ToF_Difference_For_Consider_Same_Path
                        counter = counter + 1;
                    end
                end
                NumEstPath_SpotFi = NumEstPath_SpotFi - counter;
            end
        end
        if ((NumEstPath_SpotFi > estLpaths_for_PCK(num_packet)) || (r == numel(Index_SF_search)))
            % take all the estimatin, even if the corr is less than the THR. FirstPath_AoAestimator will fix this
            matrix_Estimated_perPck(:,:,num_packet) = EstimatedMatrix(:,1:3);
            Previus_r = r;
            break
        end
    end
    
end
               
matrix_Estimated_tmp = matrix_Estimated_perPck(:,:,:); % 2D-matrix of estimation [path, [ToF, AoA, Corr]]
matrix_Estimated = matrix_Estimated_tmp(:,:);

FirstPath_AoAestimation = FirstPath_AoAestimator(matrix_Estimated,SpotFi_structure); % AoA Estimation

% This version has no moving windows. To implement it move the FirstPath_AoAestimation 
% inside the for loop and modify it appropriately, following the instruction of the Sec. 4.5
