function [OutputMatrix, MUSIC_2Dspectrum] = get_matrix_3Dinfo_SpotFi_PathAdaptation(CSI,SubCarrInd,T,c,fc,fgap,space_between_antennas,SpotFi_structure,NumEstPath_MPM)

ToF_Difference_For_Consider_Same_Path = 10;

[num_antennas, num_subcarriers] = size(CSI);
sanitizedCSI = SpotFi_Sanitization(CSI,num_subcarriers,num_antennas,SubCarrInd,T);
maxNumPaths = SpotFi_structure.max_numPaths;
[EstimatedMatrix,MUSIC_2Dspectrum] = backscatterEstimationMusic(sanitizedCSI, num_antennas, num_subcarriers, c, fc,...
    T, fgap, SubCarrInd, space_between_antennas, SpotFi_structure, ones(maxNumPaths));
%----storing high resolution temporary results----%
vector_estimatedAoA_HRSpotFi = EstimatedMatrix(1:maxNumPaths,2);
vector_estimatedToF_HRSpotFi = EstimatedMatrix(1:maxNumPaths,1);
vector_estimatedCorr_HRSpotFi = EstimatedMatrix(1:maxNumPaths,3);
MUSIC_2Dspectrum_HR = MUSIC_2Dspectrum;
%-------------------------------------------------%

vector_corrThr = EstimatedMatrix(:,3);
vector_corrThr = vector_corrThr./vector_corrThr(1);
NumEstPath_SpotFi = sum(vector_corrThr >= SpotFi_structure.selected_corrThr);

%---Modification for MPM functionality---%
ToF_pck =  vector_estimatedToF_HRSpotFi;
if NumEstPath_SpotFi > 1
    if NumEstPath_SpotFi == 2
        if abs(ToF_pck(1) - ToF_pck(2)) <= ToF_Difference_For_Consider_Same_Path
            NumEstPath_SpotFi = 1;
        end
    else
        counter = 0;
        for iiiii = 1:(NumEstPath_SpotFi - 1)
            if abs(ToF_pck(iiiii) - ToF_pck(iiiii+1)) <= ToF_Difference_For_Consider_Same_Path
                counter = counter + 1;
            end
        end
        NumEstPath_SpotFi = NumEstPath_SpotFi - counter;
    end
end
%----------------------------------------------%

NumEstPath_HRSpotFi = NumEstPath_SpotFi;
diff_NumEstPaths = NumEstPath_MPM - NumEstPath_SpotFi;
if diff_NumEstPaths == 0
    NumEstPath_firstPathEstimator = NumEstPath_MPM;
elseif diff_NumEstPaths < 0
    SpotFi_structure.L = num_subcarriers;
    SpotFi_structure.K = num_antennas;
    [EstimatedMatrix,MUSIC_2Dspectrum] = backscatterEstimationMusic(sanitizedCSI, num_antennas, num_subcarriers, c, fc,...
        T, fgap, SubCarrInd, space_between_antennas, SpotFi_structure, ones(NumEstPath_MPM));
    NumEstPath_firstPathEstimator = NumEstPath_MPM;
else
    NumIterations = 7;
    vector_L = floor(num_subcarriers:(-num_subcarriers/(2*NumIterations)):num_subcarriers/2);
    indexIteration = 2;
    while(1)
        SpotFi_structure.L = vector_L(indexIteration);
        [EstimatedMatrix,MUSIC_2Dspectrum] = backscatterEstimationMusic(sanitizedCSI, num_antennas, num_subcarriers, c, fc,...
            T, fgap, SubCarrInd, space_between_antennas, SpotFi_structure, ones(maxNumPaths));
        vector_corrThr = EstimatedMatrix(:,3);
        vector_corrThr = vector_corrThr./vector_corrThr(1);
        NumEstPath_SpotFi = sum(vector_corrThr >= SpotFi_structure.selected_corrThr);
        
        %---Modification for MPM functionality---%
        ToF_pck =  EstimatedMatrix(1:maxNumPaths,1);
        if NumEstPath_SpotFi > 1
            if NumEstPath_SpotFi == 2
                if abs(ToF_pck(1) - ToF_pck(2)) <= ToF_Difference_For_Consider_Same_Path
                    NumEstPath_SpotFi = 1;
                end
            else
                counter = 0;
                for iiiii = 1:(NumEstPath_SpotFi - 1)
                    if abs(ToF_pck(iiiii) - ToF_pck(iiiii+1)) <= ToF_Difference_For_Consider_Same_Path
                        counter = counter + 1;
                    end
                end
                NumEstPath_SpotFi = NumEstPath_SpotFi - counter;
            end
        end
        %----------------------------------------------%
        
        if NumEstPath_SpotFi >= NumEstPath_MPM || indexIteration == NumIterations
            NumEstPath_firstPathEstimator = NumEstPath_HRSpotFi;
            if NumEstPath_SpotFi == NumEstPath_HRSpotFi
                EstimatedMatrix = [vector_estimatedToF_HRSpotFi vector_estimatedAoA_HRSpotFi vector_estimatedCorr_HRSpotFi];
                MUSIC_2Dspectrum = MUSIC_2Dspectrum_HR;
            end
            break;
        end
        indexIteration = indexIteration + 1;
    end
end

OutputMatrix = EstimatedMatrix(1:NumEstPath_firstPathEstimator,1:3);

end

