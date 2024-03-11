function Testbed = get_infoTestbed(n)
Testbed = struct;

if n == 1
    Testbed.index_positions =1; 
    Testbed.totPcks_perPos = 1;
    Testbed.prefix_1 = "mat files/CSI/T1/csi_1C1_pos";
    Testbed.prefix_2 = "";
    Testbed.index_info = Testbed.index_positions;
    Testbed.MPM_PencilParameter = 5;
    Testbed.MPM_unitCircleToll = .35;
    Testbed.p_AP = [3 1.2];
    Testbed.p_UE = [.6 3.6]; 

else
    Testbed.index_positions = NaN;
    Testbed.totPcks_perPos = NaN;
    Testbed.prefix_1 = NaN;
    Testbed.prefix_2 = NaN;
    Testbed.index_info = NaN;
    Testbed.MPM_PencilParameter = NaN;
    Testbed.MPM_unitCircleToll = NaN;
    Testbed.p_AP = NaN;
    Testbed.p_UE = NaN;
end

