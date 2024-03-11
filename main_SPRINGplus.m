
% --------------------------------------------------------------------%
%| This code implements SPRING+:   
%| Please read the paper for more information about the code.   
%|                                                              
%| Authors: Stavros Eleftherakis, Maurizio Rea and Giuseppe Santaromita
%|------------------------------------------------------------------- %

%----cleaning environment----%
clear;
close all;
clc;
format
warning('off');
%----------------------------%

%----Adding paths for important folders----%
addpath('Functions');
addpath('mat files');
%------------------------------------------%

%----loading SPRING calibration mat files----%
load('delta_alpha_1');
load('delta_alpha_2');
load('delta_phi');
%--------------------------------------------%

%----Setting plots----%
plot_AoA_ECDFs = 0;
%---------------------%

%----802.11ac----%
channels_802_11ac = [7:9 11:12 16:16:32 34:2:64 68 96 100:2:128 132:2:144 149:2:155];
frequency_channels_802_11ac = 5.035e9 +((channels_802_11ac-channels_802_11ac(1))*5)*10^6;
channel = 112;
fc = frequency_channels_802_11ac(channels_802_11ac==channel);
fgap = 312.5e3; % Hz
numDataSubCarr = 234;
SubCarrInd = -122:122;
index_null_subcarriers = [-103, -75, -39, -11, -1, 0, 1, 11, 39, 75, 103];
[~,idx] = intersect(SubCarrInd,index_null_subcarriers,'stable');
SubCarrInd(idx) = [];
%----------------%

%----Fixed variables----%
c = .299792458e9; % speed of the light [m/s]
numTxAntennas = 1; % number of Tx antennas
numRxAntennas = 4; % number of Rx antennas
space_between_antennas = .026; % QTNA configuration
%-----------------------%

%----SpotFi parameters----%
SpotFi_structure = get_infoSpotFi(numDataSubCarr,numRxAntennas);
%-------------------------%

%----Setting Testbed inputs----%
selectedTestbed = 1; 
Testbed = get_infoTestbed(selectedTestbed);
%------------------------------%

estLpaths_for_PCK = NaN(max(Testbed.index_positions),Testbed.totPcks_perPos);


%-----------------------MPM application----------------------------%

for index_pos = Testbed.index_positions
    
    %----loading CSI data per position----%
    input_data(index_pos) = load(sprintf('%s%d%s',Testbed.prefix_1,Testbed.index_info(index_pos),Testbed.prefix_2));
    %-------------------------------------%    
    %----MPM algorithm----%
    estLpaths_for_PCK(index_pos,:) = MPM_delayEstimator(input_data(index_pos).H,Testbed.totPcks_perPos,numRxAntennas,Testbed.MPM_PencilParameter,Testbed.MPM_unitCircleToll);  
    %---------------------%
end



vector_FirstPath_AoAestimations = NaN(max(Testbed.index_positions));

for index_pos = Testbed.index_positions
    
     for packets=1:Testbed.totPcks_perPos
        [~,vector_FirstPath_AoAestimations(index_pos),~] = get_AoAestimations(input_data(index_pos).H,Testbed.totPcks_perPos,delta_phi,SubCarrInd,numTxAntennas,c,fc,fgap,space_between_antennas,SpotFi_structure,estLpaths_for_PCK(index_pos,:));
     end 
     
end





