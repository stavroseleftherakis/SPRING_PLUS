function estimatedAoA = FirstPath_AoAestimator(matrix_3Dinfo_SpotFi,info)

% DBSCAN Clustering and Likelihood parameters
% Need to be retuned in different scenarios
Wc = 1;
WAoA = 100;
WToF = 100;
Ws = 1;
epsilon = 0.07;
num_outliers_percent = 0.07;
num_mic_C_percent = 0.1; % For not considering clusters that contain very small number of points

singleColumn_AoA = matrix_3Dinfo_SpotFi(:,2);
singleColumn_ToF = matrix_3Dinfo_SpotFi(:,1);
singleColumn_corrThr = matrix_3Dinfo_SpotFi(:,3);
singleColumn_corrThr = singleColumn_corrThr./singleColumn_corrThr(1);
singleColumn_AoA(isnan(singleColumn_AoA)) = [];
singleColumn_ToF(isnan(singleColumn_ToF)) = [];
singleColumn_corrThr(isnan(singleColumn_corrThr)) = [];

singleColumn_ToF = singleColumn_ToF(singleColumn_corrThr>info.selected_corrThr);
singleColumn_ToF = singleColumn_ToF + abs(min(info.delayRange)*1e9);
singleColumn_ToF = singleColumn_ToF./max(singleColumn_ToF);

singleColumn_AoA = singleColumn_AoA(singleColumn_corrThr>info.selected_corrThr);
singleColumn_AoA = singleColumn_AoA./max(info.angleRange);

% Applying DBSCAN
matrix_X = [singleColumn_ToF, singleColumn_AoA];
total_points = numel(singleColumn_ToF);
num_outliers_per_cluster = max([1 floor(num_outliers_percent*total_points)]);
idx = dbscan(matrix_X,epsilon,num_outliers_per_cluster);
Num_estimated_clusters = max(idx);
C = zeros(1,Num_estimated_clusters);
Var_AoA = zeros(1,Num_estimated_clusters);
Var_ToF = zeros(1,Num_estimated_clusters);
Mean_AoA = zeros(1,Num_estimated_clusters);
Mean_ToF = zeros(1,Num_estimated_clusters);
for ii = 1:Num_estimated_clusters
    C(ii) = sum(idx == ii)./total_points;
    Var_AoA(ii) = var(singleColumn_AoA(idx==ii));
    Var_ToF(ii) = var(singleColumn_ToF(idx==ii));
    Mean_AoA(ii) = mean(singleColumn_AoA(idx==ii));
    Mean_ToF(ii) = mean(singleColumn_ToF(idx==ii));
end
if Num_estimated_clusters == 1
    estimatedAoA = Mean_AoA.*90;
else
    Likelihood = exp(Wc*C - WAoA*Var_AoA - WToF*Var_ToF - Ws*Mean_ToF);
    Likelihood(C<num_mic_C_percent) = -inf;
    estimatedAoA = Mean_AoA(Likelihood == (max(Likelihood))).*90;
end

end
