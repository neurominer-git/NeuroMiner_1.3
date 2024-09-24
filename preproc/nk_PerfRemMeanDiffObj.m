function [sY, IN] = nk_PerfRemMeanDiffObj( Y, IN )
% =========================================================================
% FORMAT [sY, IN] = nk_PerfRemMeanDiffObj( Y, IN )
% =========================================================================
% Normalizes given data (specified ind dIND) to global mean of given groups 
% (specified in sIND) by first subtracting global mean (meanY) and then 
% subtracting offsets (meanG) from target data 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2023

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y)
        % Define active indices depending on training or testing situation
        if isfield(IN,'meanY') && isfield(IN,'meanG')
           if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd{i}; else, IN.sIND =[]; end
           if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd{i}; else, IN.dIND =[]; end
        else
           if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd{i}; else, IN.sIND = []; end
           if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd{i}; else, IN.dIND = []; end
        end
        [ sY{i}, IN ] = PerfRemMeanDiffObj(Y{i}, IN ); 
    end
else
    % Define active indices depending on training or testing situation
    if isfield(IN,'meanY') && isfield(IN,'meanG')
        if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd; else, IN.sIND =[]; end
        if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd; else, IN.dIND =[]; end
    else
        if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd; else, IN.sIND = []; end
        if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd; else, IN.dIND = []; end
    end
    [ sY, IN ] = PerfRemMeanDiffObj(Y, IN );
end

% =========================================================================
function [sY, IN] = PerfRemMeanDiffObj(Y, IN)
% Defaults
global VERBOSE

if isempty(IN),eIN=true; else, eIN=false; end
sY = Y;
if isfield(IN,'func')
    mfunc = str2func(IN.func);
else
    mfunc = @nm_nanmean;
end

if eIN || ~isfield(IN,'sIND') || isempty(IN.sIND), IN.sIND = ones(size(Y,1),1); end

% Index vector / logicals matrix of subjects to apply the correction
% parameters to.
if eIN || ~isfield(IN,'dIND') || isempty(IN.dIND), IN.dIND = IN.sIND; end

% Find unque group indices and determine number of groups
if size(IN.sIND,2) > 1 % dummy coded, the case with MultiCentIntensNorm2
    MS = size(IN.sIND,2);
    sIND = zeros(size(IN.sIND,1),1);
    for i = 1:MS
        sIND(IN.sIND(:,i)) = i; 
    end
    IN.sIND = sIND; 
    MS = unique(IN.sIND); 
    % find non-zero entries in source and destination vectors
    indGS = IN.sIND~=0; indGD = IN.dIND~=0;
else
    % find non-zero entries in source and destination vectors
    indGS = IN.sIND~=0; indGD = IN.dIND~=0;
    MS = unique(IN.sIND(indGS)); 
end
nMS = numel(MS); 

if ~MS, nMS = 0; end
dIND = zeros(size(IN.dIND,1),1);
if size(IN.dIND,2) > 1
    MD = size(IN.dIND,2);
    for i = 1:MD
        dIND(IN.dIND(:,i)) = i; 
    end
    IN.dIND = dIND;
    MD = unique(IN.dIND); 
else
    MD = unique(IN.dIND(indGD)); 
end
nMD = numel(MD);  

if ~MD, nMD = 0; end

if ~isfield(IN,'meanY') || isempty(IN.meanY)  
    % Compute overall mean of the data to be adjusted in non-zero indices
    IN.meanY = mfunc(Y(indGS,:));
end
    
if ~isfield(IN,'meanG') || isempty(IN.meanG) && nMS > 0
    % find how many groups are in the source vector
    D = size(Y,2);
    % Compute group-specific means as coded in IN.sIND and store means 
    IN.meanG = zeros(nMS, D);
    for i = 1:nMS
        indGi = IN.sIND == MS(i);
        IN.meanG(i,:) = mfunc(Y(indGi,:)); 
    end
    % Store unique source groups
    IN.MS = MS;
end

% Loop through destination groups
for i=1:nMD
    
    % retrieve indices to cases of current destination group
    Di = IN.dIND == MD(i);
    
    % is the a stored group index for current destination group?
    Mi = find(IN.MS == MD(i));
    
    if Mi % yes, precomputed means exist for destination group
        if VERBOSE, fprintf('\nAdjusting destination group %g using saved means of group %g.',MD(i), MD(i)); end
        meanGi = IN.meanG(Mi,:);
       
    elseif any(MS == MD(i))
        % if not, check whether for given destination group a
        % corresponding source group exists and compute the source mean on
        % the fly
        if VERBOSE, fprintf('\nNo saved group means for destination group %g available. Computing means on the fly using the respective source group',MD(i)); end
        meanGi = mfunc(Y(IN.sIND == MD(i),:));
    else
        % if not, compute the destination-group specific mean 
        if VERBOSE, fprintf('\nNeither saved group mean or source group exists for destination group %g. Use destination group to compute means',MD(i)); end
        meanGi = mfunc(Y(IN.dIND == MD(i),:));
        if all(indGS == 0) 
            IN.meanY = nm_nanmean(Y(indGD,:)); 
        end
    end
    
    % Remove offsets between destination group and meanY.
    sY(Di,:) = Y(Di,:) - (meanGi - IN.meanY);
    
end

sY(~isfinite(sY))=0;
