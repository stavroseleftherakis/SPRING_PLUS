function sample_csi_trace_sanitized = SpotFi_Sanitization(received_CSI,N,M,SubCarrInd,T)

received_CSI = received_CSI.'; % NxM
CSI_SpotFi = reshape(received_CSI,[(M*N),1]);
sample_csi_trace_sanitized = reshape(CSI_SpotFi,N,M);
[PhsSlope, PhsCons] = removePhsSlope(sample_csi_trace_sanitized,M,SubCarrInd,N);
ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
sample_csi_trace_sanitized = sample_csi_trace_sanitized.*ToMult;
relChannel_noSlope = reshape(sample_csi_trace_sanitized, N, M, T);
sample_csi_trace_sanitized = relChannel_noSlope(:);

end