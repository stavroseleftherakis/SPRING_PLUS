function SpotFi_structure = get_infoSpotFi(N,M)
SpotFi_structure = struct;
SpotFi_structure.GridPts = [201 101 1]; % [ToF, AoA, 1]
SpotFi_structure.angleRange = [-90 90]; % range available using a ULA
SpotFi_structure.delayRange = [-200 200]*1e-9; 
SpotFi_structure.do_second_iter = 0;
SpotFi_structure.deltaRange = [0 0];
SpotFi_structure.T = 1;
SpotFi_structure.useNoise = 0;
SpotFi_structure.maxRapIters = Inf;
SpotFi_structure.generateAtot = 2;
SpotFi_structure.max_numPaths = 5;
SpotFi_structure.L = N;
SpotFi_structure.K = M;
SpotFi_structure.selected_corrThr = 0.8; % selected correlation threshold
end

