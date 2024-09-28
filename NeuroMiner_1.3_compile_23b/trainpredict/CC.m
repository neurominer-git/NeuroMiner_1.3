% =========================================================================
% FORMAT param = CC(expected, predicted)
% =========================================================================
% Compute Correlation Coefficient
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2023

function param = CC(expected, predicted)
if isempty(expected), param = []; return; end
idxnan = ~isnan(predicted) & ~isnan(expected);
param = corrcoef(expected(idxnan),predicted(idxnan));
param = param(2);
if isnan(param) 
    if numel(unique(predicted))==1
        param = 0; 
    else
        warning('Prediction algorithm returned non-finite performance measure'); 
    end
end
end