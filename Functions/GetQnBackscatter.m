
function [Pn,Ps,Qn,Qs,EigenInfo] = GetQnBackscatter(X,EigDiffCutoff, nComps)

% % by singular value
% [U,D,~] = svd(X);
% % by eigenvalue

[Utmp,D] = eig(X*X');

D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);

% % % density based clustering:
% [class,~] = DBSCAN(diag(D),2);
% SignalEndIdx = find(class==-1, 1, 'last');

minMP = 2;
useMDL = 0;
useDiffMaxVal = 0; % Default value is 1. If set to 1, it considers only those eignevalues who are above a certain threshold when compared to the maximum eigenvalue

% % % MDL criterion based
MDL = [];
lambdaTot = diag(D);
subarraySize = size(X,1);
nSegments = size(X,2);
maxMultipath = length(lambdaTot); % previously used 6 for maximum number of multipath
for k = 1:maxMultipath
    MDL(k) = -nSegments*(subarraySize-(k-1))*log(geomean(lambdaTot(k:end))/mean(lambdaTot(k:end))) + 0.5*(k-1)*(2*subarraySize-k+1)*log(nSegments);
end
% % Another attempt to take the number of multipath as minimum of MDL
[~, SignalEndIdxTmp] = min(MDL);



if useMDL
    SignalEndIdx = max(SignalEndIdxTmp-1, 1);
else
    % % % Older way of finding the SignalEndIdx based on thresholding and
    % eigenvalue difference cutoff
    VecForSignalEnd = wkeep(diag(D),5,'l'); 
    diag(D(1:floor(length(D)/2)));
    
    Criterion1 = diff(db(VecForSignalEnd))<=max(-EigDiffCutoff,min(diff(db(VecForSignalEnd))));
    Criterion3 = (VecForSignalEnd(1:end-1)/VecForSignalEnd(1)>0.03); %  Previously used 0.165 and right now using 0.065
    % % % previously used criterion
    % SignalEndIdx = find(Criterion1,1,'last');
    SignalEndIdx = find(Criterion1 & Criterion3,1,'last');
    if isempty(SignalEndIdx)
        SignalEndIdx = find(Criterion3,1,'last');
    end
end

SignalEndIdx = max(SignalEndIdx, minMP);
if ~isempty(nComps)
    SignalEndIdx = nComps; 
end



Qn = U(:,SignalEndIdx+1:end);
Pn = (Qn*Qn');
Qs = U(:,1:SignalEndIdx);
Ps = (Qs*Qs');

EigenInfo = struct;
EigenInfo.Umatrix = U;
EigenInfo.nMPs = SignalEndIdx;
EigenInfo.singularValuesdB = db(diag(D));

