function [param, model] = nk_GetParam2(Y, label, Params, ModelOnly, FeatGroups)
% =========================================================================
% FORMAT [param, model] = nk_GetParam2(Y, label, Params, ModelOnly, FeatGoups)
% =========================================================================
% Generic interface function to the algorithm-specific training modules
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2024
global TRAINFUNC SVM TIME CV CVPOS SAV DEBUG

% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label,[],'prune_single');

timevec=[];
if ~isempty(TIME) && strcmp(SVM.prog,'WBLCOX')
    % Get training index
    TrInd = CV.TrainInd{CVPOS.CV2p,CVPOS.CV2f}(CV.cvin{CVPOS.CV2p,CVPOS.CV2f}.TrainInd{CVPOS.CV1p,CVPOS.CV1f});
    if CVPOS.fFull
        TrInd = [TrInd; CV.TrainInd{CVPOS.CV2p,CVPOS.CV2f}(CV.cvin{CVPOS.CV2p,CVPOS.CV2f}.TestInd{CVPOS.CV1p,CVPOS.CV1f})];
    end
    % Extract time vector
    timevec = TIME(TrInd);
    % Recode label into event [1] vs. no event [0]
    label(label==-1)=0;
end

% Check if NaNs are in matrix and throw an error to avoid non-sense results!
sNaN = sum(~isfinite(Y(:)));
if sNaN
    if size(Y,2)<500
        writetable(table(Y), sprintf('TrainingData_Error_CV2-%g-%g_CV1-%g-%g.xlsx', CVPOS.CV2p, CVPOS.CV2f, CVPOS.CV1p, CVPOS.CV1f ))
    end
    error('\nFound %g non-finite values in training matrix!\nThis usually happens in intermediate fusion mode when some data modalities have cases with completely missing data.\nCheck your preprocessing settings and your data!', sNaN)
end

if size(Y,1) ~= size(label,1)
    error('\nTraining data matrix and labels must have the same number of observations!')
end

% Pass training matrix, labels, (and time vector) to used-defined training module
switch SVM.prog
    case 'SEQOPT'
        if  ~exist('FeatGroups','var') || isempty(FeatGroups)
            [param, model] = TRAINFUNC(Y, label, [], ModelOnly, Params );
        else
            [param, model] = TRAINFUNC(Y, label, FeatGroups, ModelOnly, Params );
        end
    case 'WBLCOX'
        [param, model] = TRAINFUNC(Y, label, timevec, ModelOnly, Params );
    otherwise
        [param, model] = TRAINFUNC(Y, label, ModelOnly, Params );
end

if ~isempty(DEBUG) && isfield(DEBUG,'eachmodel') && DEBUG.eachmodel
    filename = fullfile(pwd,sprintf('Model_%s_CV2%g-%g_CV1%g-%g.mat',SAV.matname,CVPOS.CV2p,CVPOS.CV2f,CVPOS.CV1p,CVPOS.CV1f));
    save(filename,"Y","label","model");
end
