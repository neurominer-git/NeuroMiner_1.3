function IN  = nk_PerfFeatRankObj(oY, IN)
% =========================================================================
% FORMAT function IN = nk_PerfFeatRankObj(Y, IN)
% =========================================================================
% Takes input matrix Y and ranks its features according to their relevance
% for predicting IN.curlabel. The relevance vector/matrix W can be used to
% upweight (select) or downweight (remove) respective features in the
% subsequent processing steps (see nk_PerfWActObj.m).
%
% Inputs/Outputs: 
% -------------------------------------------------------------------------
% Y                   :     M cases x N features data matrix
% IN.curlabel         :     Label vector/matrix for ranking
% IN.opt              :     Parameter optimization array
% IN.Params_desc      :     Parameter descriptions
% IN.algostr          :     Algorithm used for feature weightung
% IN.(algostr)        :     Parameter substructure for feature weighting
%                           algorithm
% IN.weightmethod     :     Upweight (=1) or downweight (=2) features
% IN.W                :     The ranking vector/matrix over the features
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2023

global VERBOSE

% Defaults
if isempty(IN),             error('Suitable input structure is missing. See the functions'' help for more information.'); end
if ~isfield(IN,'curlabel'), error('Label vector/matrix is missing! Add a ''curlabel'' variable to your input structure.'); end
if ~isfield(IN,'algostr'),  error('No feature weighting algorithm found in the input structure. Provide one according to the function help.'); end
if isfield(IN,'curglabel') && ~isempty(IN.curglabel)
    Y = oY(logical(IN.curglabel),:); L = IN.curlabel(logical(IN.curglabel)); 
else
    Y = oY; L = IN.curlabel;
end

Params_desc = []; opt =[]; 
if isfield(IN,'opt') && ~isempty(IN.opt), Params_desc = IN.Params_desc; opt = IN.opt; end

% Remove unlabeled subjects for supervised algorithms
if ~strcmp(IN.algostr,'idetect')
    indnan = isnan(L); if any(indnan), L(indnan) = []; Y(indnan,:)=[]; end
end

switch IN.algostr
    
    case 'varfeat'
        IN.W = nm_nanvar(Y);
    
    case 'idetect'
        IN.idetect.sigma = nk_ReturnParam('Sigma',Params_desc, opt); 
        IN.idetect.lambda = nk_ReturnParam('Lambda',Params_desc, opt);
        [IN.obj,IN.W] = iDetect(Y',IN.idetect); IN.W = IN.W';
        
    case 'imrelief'
        IN.imrelief.sigma = nk_ReturnParam('Sigma',Params_desc, opt); 
        IN.imrelief.lambda = nk_ReturnParam('Lambda',Params_desc, opt); 

        IN.W = IMRelief_Sigmoid_FastImple(Y', L, IN.imrelief.distance, ...
                                               IN.imrelief.sigma, ...
                                               IN.imrelief.lambda, ...
                                               IN.imrelief.maxiter, ...
                                               IN.imrelief.plotfigure, VERBOSE);
    case 'simba'
        beta = nk_ReturnParam('Beta',Params_desc, opt); 
        if ~isempty(beta)
            IN.simba.simba.extra_param.beta = beta;
        else
            IN.simba.simba.extra_param.beta = suggestBeta(Y, L);
        end

        switch IN.simba.simba.utilfunc
            case 1 % linear
                if VERBOSE, fprintf(' Simba linear'); end
            case 2 % sigmoid
                if VERBOSE, fprintf(' Simba sigmoid (beta=%g)',IN.simba.simba.extra_param.beta); end
        end

        IN.W = nk_SimbaMain(Y, L, IN.simba.simba.extra_param); IN.W = IN.W';

    case 'gflip'
        beta = nk_ReturnParam('Beta',Params_desc, opt); 
        if ~isempty(beta)
            IN.gflip.gflip.extra_param.beta = beta;
        else
            IN.gflip.gflip.extra_param.beta = suggestBeta(Y, L);
        end
        switch IN.gflip.gflip.utilfunc
            case 1
                if VERBOSE, fprintf(' G-flip zero-one'); end
            case 2
                if VERBOSE, fprintf(' G-flip linear'); end
            case 3
                if VERBOSE, fprintf(' G-flip sigmoid (beta=%g)', IN.gflip.gflip.extra_param.beta); end
        end

        [~, IN.W] = gflip(Y, L, IN.gflip.gflip.extra_param);

    case 'feast'
        if VERBOSE, fprintf(' feast'); end
        IN.FEAST.NumFeat = nk_ReturnParam('NumFeat', Params_desc, opt); 
        [~, IN.W] = nk_FEAST(Y, L, [], IN.FEAST);
        %IN.W = IN.W';
        
    case 'auc'
        % Area-under-the-Curve ooperator
        IN.W = nk_AUCFeatRank(Y, L); 

    case {'pearson','spearman'}
        % simple univariate correlation using Pearson's or Spearman's
        % correlation coefficient
        if VERBOSE, fprintf(' %s', IN.algostr); end
        if strcmp(IN.algostr,'pearson')
            type = 'pearson_fast';
        else
            type = 'spearman';
        end
        IN.W = abs(nk_CorrMat(Y,L,type));

    case 'fscore'
        % simple measure for binary classification providing a ration of
        % between-group distance to cross-group scatter. Two versions are
        % available , one based on mean/SD and the other based on
        % median/IQR. Latter one can only be accessed by adding the
        % FScoreType flag to NM.TrainParam.PREPROC{x}.ACTPARAM{y}.RANK 
        if ~isfield(IN,'FScoreType'), IN.FScoreType = 'mean'; end 
        N = []; 
        switch IN.FScoreType 
            case 'mean'
                if VERBOSE; fprintf(' F-Score [Mean/SD]'); end
                meanfun =  'nm_nanmean'; stdfun = 'nm_nanstd';
            case 'median'
                if VERBOSE; fprintf(' F-Score [Median/IQR]'); end
                meanfun = 'nm_nanmedian'; stdfun = 'iqr'; 
        end
        IN.W = nk_FScoreFeatRank(Y, L, N, meanfun, stdfun);

    case 'rgs'
        if ~isempty(opt) && ~isempty(Params_desc) 
            IN.RGS.extra_param.k = nk_ReturnParam('K',Params_desc, opt); 
            IN.RGS.extra_param.beta = nk_ReturnParam('Beta',Params_desc, opt); 
        end
        if VERBOSE, fprintf(' RGS'); end
        IN.W = RGS(Y, L, IN.RGS.extra_param)'; 

    case {'libsvm','liblin'}
        if ~isempty(opt) && ~isempty(Params_desc) 
            IN.SVM.SlackParam = nk_ReturnParam('SlackRank',Params_desc, opt); 
            IN.SVM.EpsParam = nk_ReturnParam('EpsilonRank',Params_desc, opt); 
            IN.SVM.NuParam = nk_ReturnParam('NuRank',Params_desc, opt); 
            IN.SVM.Tolerance = nk_ReturnParam('TolRank',Params_desc, opt); 
        end
        if strcmp(IN.SVM.modeflag,'regression'), L = nk_ScaleData(L,0,1); end
        if isfield(IN.SVM,'evalperf'), evalperf = IN.SVM.evalperf; else, evalperf = 1; end
        IN.W = abs(nk_SVMFeatRank(Y, L, IN.SVM, evalperf)); 

    case 'relief'
        % Currently we use the MATLAB implementation of RELIEF but we may
        % switch to the RELIEF function which does not depend on the
        % Machine Learning and Statistics toolbox
        [~,IN.W] = relieff(Y, L, IN.Relief.k);
        
    case 'anova'
        % This is a relevant filter method for multi-group settings.
        if width(L)==1
            % if zeros are detected the respective cases should be ignored
            idx_include = L~=0;
            L = nk_MakeDummyVariables(L(idx_include));
        else
            idx_include = any(L,2);
            L = L(idx_include,:);
        end
        IN.X = [ones(height(L),1) L]; 
        IN = nk_PerfANOVAObjNew(Y(idx_include,:),IN);
        IN.W = IN.F;
    
    case 'pls'
        if strcmp(IN.PLS.algostr,'spls')
            IN.PLS.cu = nk_ReturnParam('SPLS-cu',Params_desc, opt); 
            IN.PLS.cv = nk_ReturnParam('SPLS-cv',Params_desc, opt); 
        end
        [~,~,~,IN.PLS] = nk_PLS(Y, L, IN.PLS);
        IN.W = abs(IN.PLS.mpp.u);
        
    case 'extern'
        % External ranking using weight vectors generated outside of NM.
        % Use with caution: it is strongly discouraged to use weight
        % vectors that have been computed on the given NM cases outside of 
        % the cross-validation loop.
        IN.W = abs(IN.EXTERN);
		
	case 'extern2'
        % This option is not accessible in the configurator
		W1 = abs(IN.EXTERN);
		for i=1:numel(IN.algostr2)
			IN2 = IN;
			IN2.algostr = IN.algostr2{i};
			if VERBOSE, fprintf(' ...adding ranking map using %s',IN2.algostr); end
			IN2 = nk_PerfFeatRankObj(Y, IN2);
			W2 = nk_PerfScaleObj(IN2.W');
			W1 = feval(IN.operator, W1, W2');
		end
		IN.W = W1;
end

% Transpose weights if needed
if size(IN.W,1) < size(IN.W,2); IN.W = IN.W'; end

% Scale from realmin to 1 and take care of non-finite values in the weight vector
IN.W = scaledata(IN.W,[],0,1);
IN.W(IN.W==0 | isnan(IN.W)) = realmin;
IN.W(isinf(IN.W)) = 1;

% If downweighting has been selected invert the weight vector
if IN.weightmethod == 2, IN.W = 1-IN.W; end


