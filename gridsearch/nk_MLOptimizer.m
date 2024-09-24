% =========================================================================
% FORMAT GDanalysis = nk_MLOptimizer(inp, strout, id, GridAct, batchflag)
% =========================================================================
% This function performs fixed grid search optimization of machine learning
% parameters (e.g. slack, kernel params).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2024
function GDanalysis = nk_MLOptimizer(inp, strout, id, GridAct, batchflag)

global SVM RFE SAV GRD MULTI MODEFL BATCH ENSEMBLE MKLRVM CV xCV DATID CL RAND MULTILABEL PREPROC W2AVAIL OCTAVE xNM fromData CVPOS
CL = {'b*-','r*-','g*-','y*-','m*-','c*-','k*-'};
W2AVAIL = false;
   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
label           = inp.label;        % targets to predict
nclass          = inp.nclass;       % # of binary comparisons
ngroups         = inp.ngroups;      % # of groups in the analysis
lx              = inp.l;            % # of subjects
probflag        = inp.probflag;     % translate decision values to probabilities
Params          = inp.Params;
Params_desc     = inp.Params_desc;
[nl, lsel]      = nk_GetLabelDim(MULTILABEL); % Multi-label mode?
n_preml         = inp.nPreMLparams(1) - 1 ;
modalvec        = inp.ModalityVec;
if isfield(inp,'simFlag') && inp.simFlag
    [ix, jx]        = size(xCV(1).TrainInd);
else
    [ix, jx]        = size(CV(1).TrainInd);
end

% Addition 22/01/2024: Enable the use of simulated annealing for parameter
% optimization. 
if isfield(inp,'ParamOptimizationMode')
    ParamOptimizationMode = inp.ParamOptimizationMode;
else
    ParamOptimizationMode = 'bruteforce';
end

tdir = pwd; if isfield(inp,'rootdir'), tdir = inp.rootdir; end
% overwrite existing CVdatamats
ovrwrtGD = false; if isfield(inp,'ovrwrtGD'), ovrwrtGD    = inp.ovrwrtGD; end    
% overwrite existing CVresults
ovrwrtRes   = false; if isfield(inp,'ovrwrtRes'), ovrwrtRes   = inp.ovrwrtRes; end  
% save modified CVdatamat files
updGD       = false; if isfield(inp,'updGD'), updGD = inp.updGD; end

if (isfield(inp,'preprocmat') && ~isempty(inp.preprocmat))  && (~isfield(inp,'gdmat') || isempty(inp.gdmat))
    GDfl            = false;
    preprocmat      = inp.preprocmat;   % paths to preprocessed data files
elseif (isfield(inp,'gdmat') && ~isempty(inp.gdmat)) && (~isfield(inp,'preprocmat') || isempty(inp.preprocmat))
    GDfl            = true;
    gdmat           = inp.gdmat;        % paths to precomputed GD structures
elseif ( ~isfield(inp,'gdmat') || isempty(inp.gdmat)) && (~isfield(inp,'preprocmat') || isempty(inp.preprocmat))
    GDfl            = -1;
else
    GDfl            = true;
    gdmat           = inp.gdmat;
    preprocmat      = inp.preprocmat;
end

if isfield(inp,'time2event'), time2event = inp.time2event; end

if ( ~exist('batchflag','var') || isempty(batchflag)) || isempty(BATCH), batchflag = false; BATCH = false; end

detrendfl = false;

switch MODEFL
    case 'classification'
        if RAND.Decompose == 9, binmode = 0; else, binmode = 1; end    
        ngroups = numel(unique(label(~isnan(label))));
         GDanalysis.predictions = cell(lx,nclass,nl);
        if ix>1
            GDanalysis.CV2grid.predictions      = nan(lx, ix, nclass, nl);
            GDanalysis.CV2grid.CI_predictions   = nan(lx, 2, nclass);
            GDanalysis.CV2grid.mean_predictions = nan(lx, nclass, nl);
            GDanalysis.CV2grid.std_predictions  = nan(lx, nclass, nl);
            GDanalysis.CV2grid.BAC              = nan(ix, nclass, nl);
            GDanalysis.CV2grid.sens             = nan(ix, nclass, nl);
            GDanalysis.CV2grid.spec             = nan(ix, nclass, nl);
            GDanalysis.CV2grid.accuracy         = nan(ix, nclass, nl);
            GDanalysis.CV2grid.PPV              = nan(ix, nclass, nl);
            GDanalysis.CV2grid.NPV              = nan(ix, nclass, nl);
            GDanalysis.CV2grid.AUC              = nan(ix, nclass, nl);
            GDanalysis.CV2grid.DOR              = nan(ix, nclass, nl);
        end
    case 'regression'
        ngroups = 1;
        % Check whether nu- or epsilon-SVR have been selected using LIBSVM
        if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend
            detrendfl = true;
        end
        GDanalysis.predictions = cell(lx,1,nl);
        GDanalysis.CV2grid.predictions = nan(lx, ix , 1, nl);
        binmode = 0;
end
    
DISP.dimension = [];
GDanalname  = [strout '_CVresults' inp.varstr '_ID' id '.mat'];
GDanalpth   = fullfile(tdir,GDanalname);

if ~ovrwrtRes && GDfl
    if exist(GDanalpth,'file')
        if ovrwrtRes
           fprintf('\nFound CVresults-file. Overwrite it!')
        else
           fprintf('\nFound CVresults-file. Do not overwrite!')
           load(GDanalpth,'GDanalysis')
           return
        end
    end
end

%% Build parameter combination array for optimization
Ps = cell(nclass,1); nPs = zeros(nclass,1); nPdims = zeros(nclass,1); combcell = false;
for curclass = 1:nclass
    if isempty(Params) || isempty(Params{curclass})
        nPs(curclass) = 1; nPdims(curclass) = 1;
        Ps{curclass} = NaN;
    else
        Ps{curclass}  = allcomb2(Params{curclass},'matlab');
        if iscell(Ps{curclass}), combcell= true; end
        [nPs(curclass), nPdims(curclass)] = size(Ps{curclass});
    end
end

algostr = GetMLType(SVM);
DISP.algostr = algostr; DISP.figbin=[]; DISP.binwintitle = '';
DISP.Pdesc = Params_desc;

% Setup CV2 container variables:
if ~exist('GridAct','var') || isempty(GridAct), ...
        GridAct=nk_CVGridSelector(ix,jx); end

%%%%%%%%%%%%%%%%%%%%%%% SETUP GDanalysis STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%
GDanalysis.Model.ParamCombs                 = Ps;
GDanalysis.Model.ParamDesc                  = Params_desc;
GDanalysis.Model.NumParamCombs              = nPs;
GDanalysis.Model.NumParamDims               = nPdims;
GDanalysis.Model.NumPreMLParams             = n_preml + 1;
GDanalysis.Model.ModalityVec                = modalvec;
GDanalysis.GridAct                          = GridAct;
GDanalysis.RootPath                         = tdir;
GDanalysis.GDfilenames                      = cell(ix,jx);
GDanalysis.nclass                           = nclass;
GDanalysis.grid.mean_CVPerf                 = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.mean_TSPerf                 = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.mean_Err_CVTSPerf           = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.mean_Complexity             = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.mean_CVDiversity            = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.mean_TsDiversity            = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.grid.SelNodeFreq                 = nan(nPs(1),nclass,ix*jx,nl);
GDanalysis.NumModels                        = zeros(nclass,ix*jx,nl);
% Some algorithms require specific variables
switch SVM.prog
    case 'SEQOPT'
        nE = size(SVM.SEQOPT.C,2);
        GDanalysis.grid.mean_mSEQI              = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQI                = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.mean_mSEQE              = nan(nPs(1),nE,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQE                = nan(nPs(1),nE,nclass,ix*jx,nl);
        GDanalysis.grid.mean_mSEAPU             = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQAU               = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.mean_mSEQAL             = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQAL               = nan(nPs(1),nE-1,nclass,ix*jx,nl); 
        GDanalysis.grid.mean_mSEQPU             = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQPU               = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.mean_mSEQPL             = nan(nPs(1),nE-1,nclass,ix*jx,nl);
        GDanalysis.grid.sd_mSEQPL               = nan(nPs(1),nE-1,nclass,ix*jx,nl); 
        GDanalysis.grid.mean_SeqPerfGains       = nan(nPs(1),nE,nclass,ix*jx,nl);
        GDanalysis.caseprops                    = cell(lx,nclass,nl);
        GDanalysis.decvaltraj                   = cell(lx,nclass,ix,nl);
    case 'WBLCOX'
        GDanalysis.grid.mean_mCutOffPerc        = zeros(nPs(1),nclass,ix*jx,nl);
        GDanalysis.grid.mean_sdCutOffPerc       = zeros(nPs(1),nclass,ix*jx,nl);
        GDanalysis.grid.mean_mCutOffProb        = zeros(nPs(1),nclass,ix*jx,nl);
        GDanalysis.grid.mean_sdCutOffProb       = zeros(nPs(1),nclass,ix*jx,nl);
        GDanalysis.predtimes                    = cell(lx,nclass,nl);
        GDanalysis.optcutoffs                   = cell(lx,nclass,nl);
        GDanalysis.optcutoffpercs               = cell(lx,nclass,nl);
        GDanalysis.CV2cutoffs                   = cell(lx,nclass,nl);
end

% ... and multi-class must be treated separately
if MULTI.flag
    GDanalysis.multi_bestTR             = nan(ix,jx,nl);
    GDanalysis.multi_bestTS             = nan(ix,jx,nl);
    if ~MULTI.BinBind
        if isfield(GRD,'NodeSelect') && ( GRD.NodeSelect.mode == 2 || GRD.NodeSelect.mode == 3 )
            GDanalysis.multi_bestPpos   = cell(ix*jx,nl);
        else
            GDanalysis.multi_bestPpos   = zeros(ix*jx,nl);
        end
    end
    GDanalysis.multi_bestP              = cell(nclass,1);
    GDanalysis.multi_predictions        = cell(lx,nl);
    GDanalysis.multi_CV2predictions     = cell(lx,nl);
    GDanalysis.multi_probabilities      = cell(lx,ngroups,nl);
    GDanalysis.multi_CV2probabilities   = cell(lx,ngroups,nl);
    GDanalysis.grid.MultiCVPerf         = nan(nPs(1),ix*jx,nl);
    GDanalysis.grid.MultiTSPerf         = nan(nPs(1),ix*jx,nl);
    GDanalysis.grid.MultiERR_CVTSPerf   = nan(nPs(1),ix*jx,nl);
    GDanalysis.grid.MultiCVDiversity    = nan(nPs(1),ix*jx,nl);
    GDanalysis.grid.MultiTsDiversity    = nan(nPs(1),ix*jx,nl);
    GDanalysis.grid.MultiSelNodeFreq    = nan(nPs(1),ix*jx,nl);
    DISP.figmulti = [];    
end

GDanalysis.bestTR                           = cell(nclass,1);
GDanalysis.bestTS                           = cell(nclass,1);
GDanalysis.bestP                            = cell(nclass,1);
GDanalysis.bestPpos                         = cell(nclass,1);
GDanalysis.bestComplexity                   = cell(nclass,1);
GDanalysis.bestError                        = cell(nclass,1);

for h=1:nclass
    
    GDanalysis.bestTR{h}                    = zeros(ix,jx,nl);
    GDanalysis.bestTS{h}                    = zeros(ix,jx,nl);
    GDanalysis.bestComplexity{h}            = zeros(ix,jx,nl);
    GDanalysis.bestError{h}                 = zeros(ix,jx,nl);
    
    if (isfield(GRD,'NodeSelect') &&  GRD.NodeSelect.mode ~= 1) || combcell
        GDanalysis.bestP{h}                 = cell(ix*jx,nl);
        if MULTI.flag, GDanalysis.multi_bestP{h} = cell(ix*jx,nl); end
        GDanalysis.bestPpos{h}              = cell(ix*jx,nl);
    else
        GDanalysis.bestP{h}                 = zeros(ix*jx,nPdims(h),nl);
        if MULTI.flag, GDanalysis.multi_bestP{h} = zeros(ix*jx,nPdims(h),nl); end
        GDanalysis.bestPpos{h}              = zeros(ix*jx,nl);
    end
end
[~, ~, ~, ~, act] = nk_ReturnEvalOperator(SVM.GridParam);

if ~batchflag && RFE.dispres
   DISP.binwintitle = sprintf('NM Optimization Viewer => Analysis [#%g]: %s', inp.curanal, inp.P.analysis_id);
end

% Parameter flag structure for preprocessing
paramfl         = struct('use_exist',true,'found', false, 'write', true);

ol = 0; ll = 1; GridUsed = false(size(GridAct)); 

% If we are operating on the fly, we potentially need to run a number of
% preprocessing operations at the single-subject level
if GDfl == -1
    
    % Scale the labels
    label = nk_LabelTransform(PREPROC,MODEFL,label);
    inp.labels = label;
    
    % Perform spatial ops on the data (imaging only)
    Y = nk_PerfSpatFilt( inp.Y, PREPROC, inp.P.X );

    % ... the calibration data (Clara), too
    if isfield(inp, 'C') && inp.C{1,1}.calibflag
        CYfile = inp.C{1,1}.Y; 
        load(CYfile, 'CY');
        inp.C{1,1}.Y = CY;
        inp.C{1,1}.Y = nk_PerfSpatFilt(inp.C{1,1}.Y, PREPROC, inp.P.X);
    end
    % ... and alternative training data (Nikos; not accessible via the NM
    % menu configuration)
    if isfield(inp,'AltY')
        inp.AltY = nk_PerfSpatFilt(inp.AltY, PREPROC, inp.P.X);
    end

    % Perform required spatial ops on mask image
    inp = nk_SmoothMask( PREPROC, inp.P );

    % Set template flag
    if isfield(PREPROC,'TEMPLPROC') && ~isempty(PREPROC.TEMPLPROC) && PREPROC.TEMPLPROC
        paramfl.templateflag = true;
    end
else
     % Scale the labels
    label = nk_LabelTransform(PREPROC,MODEFL,label);
    inp.labels = label;
end

for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds

        if ~GridAct(f,d), ll = ll +1; continue, end
        DISP.f = f; DISP.d = d;
        if isfield(inp,'simFlag') && inp.simFlag
            TsInd = xCV(1).TestInd{f,d};
        else
            TsInd = CV(1).TestInd{f,d};
        end
        
        GDxfl = false;
        inp.ll=ll;
        cvstr = ['_oCV' num2str(f) '.' num2str(d) ];
        oCVfile = [strout '_CVdatamat' cvstr inp.varstr '_ID' id '.mat'];
        oCVpath = fullfile(tdir,oCVfile);
        divstr = repmat('-',1,length(oCVfile)); fprintf('\n%s',divstr)
        CVPOS.CV2p = f;
        CVPOS.CV2f = d;

        %%%%%%%%%%%%%%%%%%%%%%% FILE CONTROL SECTION %%%%%%%%%%%%%%%%%%%%%%
        switch GDfl 
            % ---------------- COMPUTE FROM SCRATCH -----------------------
            case -1 % No PreprocData / CVdatamat specified
                if fromData
                    GD=[];
                else
                    GD=[]; if ~inp.ovrwrtGD, GD = nk_CheckLoadFile(oCVpath, 'CVdatamat', f, d, ovrwrtGD, nclass); end
                end
                if isempty(GD)
                   
                    fprintf('\nComputing PreprocData on the fly.')
                    inp.f = f; inp.d = d; inp.nclass = nclass;
                    if strcmp(MODEFL,'classification') && size(label,2)>1
                        if isfield(inp,'simFlag') && inp.simFlag
                            TCV = xCV;
                        else
                            TCV = CV;
                        end
                        for cl = 1:nl
                            curlabel = lsel(cl);
                            if isfield(inp,'simFlag') && inp.simFlag
                                xCV = TCV(curlabel);
                            else
                                CV = TCV(curlabel);
                            end
                            inp.curlabel=curlabel;
                            if inp.stacking
                                mapY(curlabel) = nk_PerfPreprocessMeta(inp, label(:,curlabel), paramfl);
                            else
                                mapY(curlabel) = nk_PerfPreprocess(Y, inp, label(:,curlabel), paramfl);
                            end
                        end
                        if isfield(inp,'simFlag') && inp.simFlag
                            xCV = TCV;
                        else
                            CV = TCV;
                        end
                    else
                        inp.curlabel = 1;
                        if inp.stacking
                            mapY = nk_PerfPreprocessMeta(inp, label(:,lsel), paramfl);
                        else
                            mapY = nk_PerfPreprocess(Y, inp, label(:,lsel), paramfl);
                        end
                        
                    end
                    ol=ol+1; GridUsed(f,d) = true;
                    
                else
                    GDanalysis.RootPath = fileparts(oCVpath);
                    GDxfl = true;
                end
                
                
            % -------------- USE PRECOMPUTED PREPROCMAT -------------------
            case 0 % PreprocData specified, but no CVdatamats
                GD = nk_CheckLoadFile(oCVpath, 'CVdatamat', f, d, ovrwrtGD, nclass);
                
                if isempty(GD)
                    ppath = preprocmat{f,d};
                    mapY = nk_CheckLoadFile(ppath, 'PreprocData', f, d, [],nclass);
                    if isempty(mapY)
                        fprintf('\n')
                        warning(['No valid preprocessed data detected for CV2 partition. Continue ' ...
                        '[' num2str(f) ', ' num2str(d) ']!']); ll=ll+1; continue
                    end
                else
                    GDanalysis.RootPath = fileparts(oCVpath);
                    GDxfl = true;
                end

                ol=ol+1; GridUsed(f,d) = true;

            case 1 % use already computed CVdatamats

                gdpath = gdmat{inp.P.curmodal}{f,d}(1,:);
                if isempty(gdpath) || ~exist(gdpath,'file')
                    fprintf('\n')
                    warning(['No valid CV2datamat detected for CV2 partition. Continue ' ...
                            '[' num2str(f) ', ' num2str(d) ']!']); ll=ll+1;
                    continue; 
                end
                p = fileparts(deblank(gdpath));
                cvstr = ['_oCV' num2str(f) '.' num2str(d)];
                oCVpath = fullfile(p,[strout '_CVdatamat' cvstr inp.varstr '_ID' id '.mat']);

                [GD, GDfound] = nk_CheckLoadFile(gdpath, 'CVdatamat', f, d, ovrwrtGD, nclass);
                if isempty(GD) && ~GDfound 
                    fprintf('\n')
                    warning(['No valid CVdatamat detected for CV2 partition ' ...
                        '[' num2str(f) ', ' num2str(d) ']!']); ll=ll+1;
                    continue
                else
                    GDanalysis.RootPath = fileparts(oCVpath);
                    GDxfl = true;
                end

                ol = ol+1; GridUsed(f,d) = true;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%% INNER LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (~GDfl || GDfl == -1) && ~GDxfl

            % Check for multiple non-concatenated variates
            try
                nvar = size(mapY.Tr,3);
            catch
                nvar = size(mapY(1).Tr,3);
            end
            
            if strcmp(SVM.prog,'MKLRVM') && nvar > 1 && ll==1
                SVM.kernel.kernstr = repmat({SVM.kernel.kernstr},1,nvar);
                MKLRVM.standardize_flag = repmat({MKLRVM.standardize_flag},1,nvar);
            end

            % %%%%%%%%%%%%%%%%%%%%% PREPARATIONS %%%%%%%%%%%%%%%%%%%%%%
            % CV1 test data performance measures
            GD.TR       = zeros(nPs(1),nclass,nl);
            GD.sTR     = zeros(nPs(1),nclass,nl);

            % CV2 test data performance measures
            GD.TS       = zeros(nPs(1),nclass,nl); 
            GD.mTS      = zeros(nPs(1),nclass,nl); % mean ?
            GD.sTS      = zeros(nPs(1),nclass,nl); % sd ?

            % Generalization error between CV1 and CV2 test data
            GD.ERR      = zeros(nPs(1),nclass,nl);

            % Final binary classifier / predictor results on CV2 test
            % data
            GD.BinPred  = cell(nPs(1),nl);

            % Multi-group classification measures
            if MULTI.flag
               GD.MultiTR       =  zeros(nPs(1),nl); % performance on CV1 test data
               GD.MultiTS       =  zeros(nPs(1),nl); % performance on CV2 test data
               GD.MultiERR      = zeros(nPs(1),nl); % generalization error
               GD.MultiCV1TrPred= cell(nPs(1),nl); % CV1 traindata predictions
               GD.MultiCV1CVPred= cell(nPs(1),nl); % CV1 test data predictions
               GD.MultiCV1TrProb= cell(nPs(1),ngroups,nl); % CV1 traindata predictions
               GD.MultiCV1CVProb= cell(nPs(1),ngroups,nl); % CV1 test data predictions
               GD.MultiPred     = cell(nPs(1),nl); % CV2 test data predictions
               GD.MultiM_DivT   = zeros(nPs(1),nl);
               GD.MultiSD_DivT  = zeros(nPs(1),nl);             
               GD.MultiM_DivV   = zeros(nPs(1),nl);
               GD.MultiSD_DivV  = zeros(nPs(1),nl);
               GD.MultiCV2Div   = zeros(nPs(1),nl);
               GD.MultiCV2DivDec= zeros(nPs(1),nl);
            end

            % Mean model complexity across CV1 partitions
            GD.C        = zeros(nPs(1),nclass,nl);

            % Diversity measures for CV1 and CV2 test data
            GD.M_DivT   = zeros(nPs(1),nclass,nl);
            GD.SD_DivT  = zeros(nPs(1),nclass,nl);             
            GD.M_DivV   = zeros(nPs(1),nclass,nl);
            GD.SD_DivV  = zeros(nPs(1),nclass,nl);
            GD.CV2Div   = zeros(nPs(1),nclass,nl);
            GD.CV2DivDec= zeros(nPs(1),nclass,nl);

            % Models params
            MD          = cell(nPs(1),nl);  % models
            GD.FEAT     = cell(nPs(1),nl);  % selected features for model in MD
            GD.VI       = cell(nPs(1),nl);
            GD.Weights  = cell(nPs(1),nl);  % weights for base learners' predictions

            % Decision values / Probabilities of CV1 training & test data and CV2 test data
            GD.DT       = cell(nPs(1),nl);  % CV1 training data
            GD.DV       = cell(nPs(1),nl);  % CV1 test data          
            GD.DS       = cell(nPs(1),nl);  % CV2 test data
            
            % For sequence optimizer only
            switch SVM.prog
               case 'SEQOPT'
                   GD.mSEQI = cell(nPs(1), nclass, nl);
                   GD.sdSEQI = cell(nPs(1),nclass, nl);
                   GD.mSEQE = cell(nPs(1), nclass, nl);  
                   GD.sdSEQE = cell(nPs(1), nclass, nl); 
                   GD.mSEQAbsThrU = cell(nPs(1), nclass, nl);
                   GD.sdSEQAbsThrU = cell(nPs(1), nclass, nl);
                   GD.mSEQAbsThrL = cell(nPs(1), nclass, nl);
                   GD.sdSEQAbsThrL = cell(nPs(1), nclass, nl);
                   GD.mSEQPercThrU = cell(nPs(1), nclass, nl);
                   GD.sdSEQPercThrU = cell(nPs(1), nclass, nl);
                   GD.mSEQPercThrL = cell(nPs(1), nclass, nl);
                   GD.sdSEQPercThrL = cell(nPs(1), nclass, nl);
                   GD.CasePropagations = cell(nPs(1), nclass, nl);
                   GD.DecValTraj = cell(nPs(1), nclass, nl);
                case 'WBLCOX'
                   GD.mCutOffPerc = zeros(nPs(1), nclass, nl);
                   GD.sdCutOffPerc = zeros(nPs(1),nclass, nl);
                   GD.mCutOffProb = zeros(nPs(1), nclass, nl);  
                   GD.sdCutOffProb = zeros(nPs(1), nclass, nl);
                   GD.CV2Cutoffs = cell(nPs(1), nclass, nl);
                   GD.CV1predictedtimes = cell(nPs(1), nclass, nl);
                   GD.CV2predictedtimes = cell(nPs(1), nclass, nl);
            end

            if detrendfl, GD.Detrend = cell(nPs(1),nl); end
            if isfield(RFE.Wrapper,'optflag') && RFE.Wrapper.optflag == 1, RFE.Wrapper.flag = 0; end
            
            %%%%%%%%%%%%%%%% PARAMETER OPTIMIZATION %%%%%%%%%%%%%%%%
            % Addition 22/01/2024: Enable the use of simulated annealing for parameter optimization. 
            switch ParamOptimizationMode
                case 'bruteforce'
                    [ GD, MD ] = nk_MLOptimizer_ParamCycler(GD, MD, DISP, Ps, Params_desc, mapY, algostr, f, d, n_preml, nclass, ngroups, batchflag, [], combcell);
                case 'simanneal'
                    [ GD, MD ] = nk_MLOptimizer_ParamAnnealer(GD, MD, DISP, Ps, Params_desc, mapY, algostr, f, d, n_preml, nclass, batchflag, [], combcell);
            end
            
            %%%%%%%%%%%%%%%% MODEL SELECTION LOGIC %%%%%%%%%%%%%%%%%
            [GD, MultiBinBind] = nk_ModelNodeSelector(GD, MD, label, f, d, nclass, Ps, Params_desc, combcell, act);
            
            %%%%%%%%%% WRAPPER-BASED LEARNING AT OPTIMA %%%%%%%%%%%% 
            if isfield(RFE.Wrapper,'optflag') && RFE.Wrapper.optflag == 1
                RFE.Wrapper.flag = 1;
                switch ParamOptimizationMode
                    case 'bruteforce'
                        [ GD, MD ] = nk_MLOptimizer_ParamCycler(GD, MD, DISP, Ps, Params_desc, mapY, algostr, f, d, n_preml, nclass, ngroups, batchflag, GD.BinaryGridSelection, combcell);
                    case 'simanneal'
                        [ GD, MD ] = nk_MLOptimizer_ParamAnnealer(GD, MD, DISP, Ps, Params_desc, mapY, algostr, f, d, n_preml, nclass, batchflag, [], combcell);
                end
                [ GD, MultiBinBind ] = nk_ModelNodeSelector(GD, MD, label, f, d, nclass, Ps, Params_desc, combcell, act);
            end
            
            if inp.stacking
                GD.nM_cnt = mapY.nM_cnt;
            end
            
            % Set saving flag to store GD on hard disk
            saveGDflag = true;
            

        else %Pre-existing data will be used

            tGD = GD;

            if ~exist('MD','var'), MD = cell(nPs(1),nl); end

            % if no mapY has been transmitted ask user
            if ~exist('mapY','var')
                if probflag 
                    fprintf('\nPreprocessed data needed for computation of probabilities!')
                    preprocmat = nk_GenPreprocMaster(DATID, CV);
                    ppath = preprocmat{dim_index,f,d};
                    mapY = nk_CheckLoadFile(ppath, 'PreprocData', dimension, f, d, [] ,nclass);
                else
                    mapY =[]; mapYi = [];
                end
            end
            %%%%%%%%%%%%%%%% SELECT PARAMS AT OPTIMUM %%%%%%%%%%%%%%%%%
            [tGD, MultiBinBind] = nk_ModelNodeSelector(tGD, MD, label, f, d, nclass, Ps, Params_desc, combcell, act);

             % Set saving flag to store GD on hard disk only if tGD ~= GD:
            if isequaln(tGD,GD) 
                saveGDflag = false; 
            else
                saveGDflag = true;
            end
            GD = tGD; clear tGD;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%% Save CVDATAMAT %%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(inp,'simFlag') && inp.simFlag
            saveGDflag = 0;
        end

        [~,oCVnam] = fileparts(oCVpath); 
        if saveGDflag && updGD
            operm = d; ofold = f; 
            if ~batchflag
                if ~exist(oCVpath,'file')
                    savflag = 1;
                else
                    if ~exist('keepflag','var')
                        savflag = nk_input(['Update ' oCVnam],0,'yes|no',[1,0],1);
                        keepflag = nk_input('Apply to all other cases',0,'yes|no',[1,0],1);
                    else
                        if ~keepflag, savflag = nk_input(['Update ' oCVnam],0,'yes|no',[1,0],1); end
                    end
                    if savflag, fprintf('\nUpdating: %s', oCVnam), end
                end
            else
                savflag = 1;
            end
            if savflag
                fprintf('\nSaving: %s.', oCVpath)
                if SAV.savemodel
                    if OCTAVE
                        save(oCVpath,'GD','MD','Ps','Params_desc','operm','ofold');
                    else
                        save(oCVpath,'-v7.3','GD','MD','Ps','Params_desc','operm','ofold');
                    end
                else
                    if OCTAVE
                        save(oCVpath,'GD','Ps','Params_desc','operm','ofold');
                    else
                        save(oCVpath,'-v7.3','GD','Ps','Params_desc','operm','ofold');
                    end
                end
                
            end
        else
            fprintf('\nUpdate of %s skipped!',oCVnam)
        end

        if batchflag == 0 || batchflag == 2
            GDanalysis.GDfilenames{f,d} = [ oCVnam '.mat' ]; 
        end
        
        if GDfl || batchflag == 0 || batchflag == 2

            %%%%%%%%%%%%%%%%%%%% CONSTRUCT GDanalysis %%%%%%%%%%%%%%%%%%%%%                

            % BINARY CLASSIFICATION or REGRESSION PERFORMANCE
            % =============================================================
            % CV1 test data performance
            GDanalysis.grid.mean_CVPerf(:,:,ll,:)         = GD.TR;
            % CV2 test data performance
            GDanalysis.grid.mean_TSPerf(:,:,ll,:)         = GD.TS;
            % Generalization error
            GDanalysis.grid.mean_Err_CVTSPerf(:,:,ll,:)   = GD.ERR;
            % Model complexity
            GDanalysis.grid.mean_Complexity(:,:,ll,:)     = GD.C;
            % Ensemble diversity CV1 test data
            GDanalysis.grid.mean_CVDiversity(:,:,ll,:)    = GD.M_DivV;
            % Ensemble diversity CV2 test data
            GDanalysis.grid.mean_TsDiversity(:,:,ll,:)    = GD.CV2Div;
            % Specifically treat the sequence optimizer algorithm
            switch SVM.prog
                case 'SEQOPT'
                    for curclass = 1 : nclass
                        % Mean sequence gain
                        GDanalysis.grid.mean_mSEQI(:,:,curclass, ll,:)   = cell2matpadnan(GD.mSEQI(:,curclass));
                        % SD sequence gain
                        GDanalysis.grid.sd_mSEQI(:,:,curclass, ll,:)     = cell2matpadnan(GD.sdSEQI(:,curclass));
                        % Mean examination frequencies
                        GDanalysis.grid.mean_mSEQE(:,:,curclass,ll,:)    = cell2matpadnan(GD.mSEQE(:,curclass));
                        % SD examination frequencies
                        GDanalysis.grid.sd_mSEQE(:,:,curclass,ll,:)      = cell2matpadnan(GD.sdSEQE(:,curclass));
                        % Mean upper percentile for case propagation
                        GDanalysis.grid.mean_mSEQAU(:,:,curclass,ll,:)   = cell2matpadnan(GD.mSEQAbsThrU(:,curclass));
                        % SD upper percentile for case propagation
                        GDanalysis.grid.sd_mSEQAU(:,:,curclass,ll,:)     = cell2matpadnan(GD.sdSEQAbsThrU(:,curclass));
                        % Mean lower percentile for case propagation
                        GDanalysis.grid.mean_mSEQAL(:,:,curclass,ll,:)   = cell2matpadnan(GD.mSEQAbsThrL(:,curclass));
                        % SD lower percentile for case propagation
                        GDanalysis.grid.sd_mSEQAL(:,:,curclass,ll,:)   = cell2matpadnan(GD.sdSEQAbsThrL(:,curclass));
                        % Mean upper percentile for case propagation
                        GDanalysis.grid.mean_mSEQPU(:,:,curclass,ll,:)   = cell2matpadnan(GD.mSEQPercThrU(:,curclass));
                        % SD upper percentile for case propagation
                        GDanalysis.grid.sd_mSEQPU(:,:,curclass,ll,:)     = cell2matpadnan(GD.sdSEQPercThrU(:,curclass));
                        % Mean lower percentile for case propagation
                        GDanalysis.grid.mean_mSEQPL(:,:,curclass,ll,:)   = cell2matpadnan(GD.mSEQPercThrL(:,curclass));
                        % SD lower percentile for case propagation
                        GDanalysis.grid.sd_mSEQPL(:,:,curclass,ll,:)   = cell2matpadnan(GD.sdSEQPercThrL(:,curclass));
                        
                        GDanalysis.grid.mean_SeqPerfGains(:,:,curclass,ll,:) = cell2matpadnan(cellfun(@nm_nanmean, GD.SeqPerfIncreases(:,curclass), 'UniformOutput', false), [], nE);
                    end
                case 'WBLCOX'
                    for curclass = 1 : nclass
                        GDanalysis.grid.mean_mCutOffPerc(:,curclass,ll) = GD.mCutOffPerc(:,curclass);
                        GDanalysis.grid.mean_sdCutOffPerc(:,curclass,ll) = GD.sdCutOffPerc(:,curclass);
                        GDanalysis.grid.mean_mCutOffProb(:,curclass,ll) = GD.mCutOffProb(:,curclass);
                        GDanalysis.grid.mean_sdCutOffProb(:,curclass,ll) = GD.sdCutOffProb(:,curclass);
                    end
            end
            
            if MULTI.flag
               GDanalysis.grid.MultiCVDiversity(:,ll,:)   = GD.MultiM_DivV; 
               GDanalysis.grid.MultiTsDiversity(:,ll,:)   = GD.MultiCV2Div; 
            end

            for curclass=1:nclass
                
                if MULTILABEL.flag
                    if isfield(MULTILABEL,'sel')
                        nl = numel(MULTILABEL.sel);
                    else
                        nl = MULTILABEL.dim;
                    end
                else
                    nl = 1;
                end
                for curlabel=1:nl

                    GDanalysis.grid.SelNodeFreq(:,curclass,ll,curlabel) = GD.BinaryGridSelection{curclass}{curlabel}.SelNodes;

                    % Best performance measures
                    if (isfield(GRD,'NodeSelect') && ( GRD.NodeSelect.mode ~= 1)) || combcell
                        GDanalysis.bestP{curclass}{ll,curlabel}      = GD.BinaryGridSelection{curclass}{curlabel}.bestP;
                        if MULTI.flag
                            GDanalysis.multi_bestP{curclass}{ll,curlabel} = GD.MultiGroupGridSelection{curlabel}.bestP{curclass};
                        end
                        GDanalysis.bestPpos{curclass}{ll,curlabel}   = GD.BinaryGridSelection{curclass}{curlabel}.Npos;
                        GDanalysis.bestTR{curclass}(f,d,curlabel)    = nm_nanmean(GD.BinaryGridSelection{curclass}{curlabel}.bestacc);
                        GDanalysis.bestTS{curclass}(f,d,curlabel)    = nm_nanmean(GD.BinaryGridSelection{curclass}{curlabel}.besttestparam);
                        GDanalysis.bestComplexity{curclass}(f,d,curlabel) = nm_nanmean(GD.BinaryGridSelection{curclass}{curlabel}.bestcomplexity);
                        GDanalysis.bestError{curclass}(f,d,curlabel) = nm_nanmean(GD.BinaryGridSelection{curclass}{curlabel}.besterr);

                        EnsDat=[];
                        for zu = 1:GD.BinaryGridSelection{curclass}{curlabel}.Nodes
                            EnsDat = [EnsDat nk_cellcat(GD.BinaryGridSelection{curclass}{curlabel}.bestpred{zu},[],2)];
                        end
                    else
                        GDanalysis.bestP{curclass}(ll,:,curlabel) = GD.BinaryGridSelection{curclass}{curlabel}.bestP(1,:);
                        if MULTI.flag
                            GDanalysis.multi_bestP{curclass}(ll,:,curlabel) = GD.MultiGroupGridSelection{curlabel}.bestP{curclass};
                        end
                        GDanalysis.bestPpos{curclass}(ll,curlabel)   = GD.BinaryGridSelection{curclass}{curlabel}.Npos(1);
                        GDanalysis.bestTR{curclass}(f,d,curlabel)    = GD.BinaryGridSelection{curclass}{curlabel}.bestacc(1);
                        GDanalysis.bestTS{curclass}(f,d,curlabel)    = GD.BinaryGridSelection{curclass}{curlabel}.besttestparam(1);
                        GDanalysis.bestComplexity{curclass}(f,d,curlabel) = GD.BinaryGridSelection{curclass}{curlabel}.bestcomplexity(1);
                        GDanalysis.bestError{curclass}(f,d,curlabel) = GD.BinaryGridSelection{curclass}{curlabel}.besterr(1);
                        % Concatenate cell array of decisions into ensemble
                        % decision matrix
                        EnsDat = nk_cellcat(GD.BinaryGridSelection{curclass}{curlabel}.bestpred{1},[],2);
                    end
                                            
                    switch binmode
                        case 0 % MULTI-GROUP or regression analysis
                            TsI = TsInd;
                        case 1 % BINARY-GROUP analysis
                            if strcmp(MODEFL,'classification')
                                if isfield(inp,'simFlag') && inp.simFlag
                                    TsI = TsInd(xCV(curlabel).classnew{f,d}{curclass}.ind); 
                                else
                                    TsI = TsInd(CV(curlabel).classnew{f,d}{curclass}.ind); 
                                end
                            else
                                if isfield(inp,'simFlag') && inp.simFlag
                                    TsI = TsInd(xCV.classnew{f,d}{curclass}.ind); 
                                else
                                    TsI = TsInd(CV.classnew{f,d}{curclass}.ind); 
                                end
                            end
                    end

                    switch MODEFL
                        case 'classification'
                            if RAND.Decompose ~=9
                                if isfield(inp,'simFlag') && inp.simFlag
                                    binInd = xCV(curlabel).classnew{f,d}{curclass}.ind;
                                else
                                   binInd = CV(curlabel).classnew{f,d}{curclass}.ind;
                                end
                            else
                                binInd = 1:size(EnsDat,1);
                            end
                        case 'regression'
                            binInd = 1:size(EnsDat,1);
                    end
                    
                    if ix>1
                        if strcmp(MODEFL,'classification')
                            % Majority voting at the CV2 level
                            GDanalysis.CV2grid.predictions(TsI, f, curclass, curlabel) = nm_nansum(sign(EnsDat(binInd,:))>0,2)./sum(~isnan(EnsDat(binInd,:)),2);
                        else
                            % Just the median for regression
                            GDanalysis.CV2grid.predictions(TsI, f, curclass, curlabel) = nm_nanmedian(EnsDat(binInd,:),2);
                        end
                    end
                    GDanalysis.NumModels(curclass,ll,curlabel) = GDanalysis.NumModels(curclass,ll,curlabel) + size(EnsDat,2); 
                    
                    % Choose Metric (hard or soft predictions):
                    % Metric = 1 => use hard labels for aggregation
                    % Metric = 2 => use decision scores/probability for aggregation
                    if RFE.CV2Class.EnsembleStrategy.Metric == 1 && ...
                            ~strcmp(MODEFL,'regression') && RAND.Decompose ~=9        
                            EnsDat = sign(EnsDat); 
                    end
                    % Check aggregation level:
                    % 0 = Mean of CV1 ensemble decision (grand mean approach)
                    % 1 = Concatenate CV1 base learners' decision into big ensemble!
                    if size(EnsDat,2)>1 && ~RFE.CV2Class.EnsembleStrategy.AggregationLevel 
                        EnsDat = nm_nanmedian(EnsDat,2);
                    end
                    
                    % Concatenate (averaged) CV1 ensemble decisions along the 
                    % column dimension for each hold-out CV2 test sample
                    GDanalysis.predictions(TsI, curclass, curlabel) = ...
                        cellmat_mergecols(GDanalysis.predictions(TsI, curclass,curlabel), ...
                                            num2cell(EnsDat(binInd,:),2));

                    % Concatenate node propagation indices in case the
                    % sequence optimization algorithm has been run on 
                    if any(strcmp(SVM.prog,{'SEQOPT','WBLCOX'}))
                        fSelNodes = find(GD.BinaryGridSelection{curclass}{curlabel}.SelNodes);
                        cv2lx = size(GD.BinaryGridSelection{curclass}{curlabel}.bestpred{1}{1,1},1);        
                        switch SVM.prog
                            case 'SEQOPT'
                                if numel(fSelNodes)>1
                                    CaseProps = zeros(cv2lx, numel(fSelNodes));
                                    DecValTraj = nan([cv2lx size(SVM.SEQOPT.C,2) numel(fSelNodes)]);
                                    for cp=1:numel(fSelNodes)
                                         CaseProps(:,cp) = nm_nanmedian(GD.CasePropagations{fSelNodes(cp)},2);
                                         DecValTraj_cp = nm_nanmedian(GD.DecValTraj{fSelNodes(cp)},3);
                                         lenTraj_cp = size(DecValTraj_cp,2);
                                         DecValTraj(:,1:lenTraj_cp,cp) = DecValTraj_cp; 
                                    end
                                    mCaseProps = nm_nanmedian(CaseProps,2);
                                    mDecValTraj = nm_nanmedian(DecValTraj,3);
                                else
                                    mCaseProps = GD.CasePropagations{fSelNodes};
                                    mDecValTraj = nm_nanmedian(GD.DecValTraj{fSelNodes},3);
                                end
                                GDanalysis.caseprops(TsI, curclass, curlabel) = cellmat_mergecols(GDanalysis.caseprops(TsI, curclass,curlabel), num2cell(mCaseProps,2));
                                GDanalysis.decvaltraj(TsI, curclass, f, curlabel) = cellmat_mergecols(GDanalysis.decvaltraj(TsI, curclass, f, curlabel), num2cell(mDecValTraj,2));
                            case 'WBLCOX'
                                if numel(fSelNodes)>1
                                    PredTimes = []; CV2Cutoffs = [];
                                    OptCutOffs = zeros(1, numel(fSelNodes));
                                    OptCutOffPercs = zeros(1, numel(fSelNodes));
                                    for cp = 1:numel(fSelNodes)
                                        PredTimes = [PredTimes nk_cellcat(GD.CV2predictedtimes{fSelNodes(cp)},[],2)];
                                        CV2Cutoffs = [CV2Cutoffs GD.CV2Cutoffs{fSelNodes(cp)}(:)'];
                                        OptCutOffs(cp) = GD.mCutOffProb(fSelNodes(cp));
                                        OptCutOffPercs(cp) = GD.mCutOffPerc(fSelNodes(cp));
                                    end
                                    mPredTimes = nm_nanmedian(PredTimes,2);
                                    mCV2Cutoffs = nm_nanmedian(CV2Cutoffs);
                                    mOptCutoffs = nm_nanmedian(OptCutOffs);
                                    mOptCutoffPercs = nm_nanmedian(OptCutOffPercs);
                                else
                                    mPredTimes = nk_cellcat(GD.CV2predictedtimes{fSelNodes(1)},[],2);
                                    mCV2Cutoffs = nm_nanmedian(GD.CV2Cutoffs{fSelNodes(1)}(:));
                                    mOptCutoffs = GD.mCutOffProb(fSelNodes(1));
                                    mOptCutoffPercs = GD.mCutOffPerc(fSelNodes(1));
                                end
                                GDanalysis.CV2cutoffs(TsI, curclass, curlabel) = cellmat_mergecols( GDanalysis.CV2cutoffs(TsI, curclass, curlabel), num2cell(repmat(mCV2Cutoffs, numel(TsI),1)));  
                                GDanalysis.predtimes(TsI, curclass, curlabel) = cellmat_mergecols( GDanalysis.predtimes(TsI, curclass, curlabel), num2cell(mPredTimes,2));  
                                GDanalysis.optcutoffs(TsI, curclass, curlabel) = cellmat_mergecols( GDanalysis.optcutoffs(TsI, curclass, curlabel), num2cell(repmat(mOptCutoffs, numel(TsI),1)));  
                                GDanalysis.optcutoffpercs(TsI, curclass, curlabel) = cellmat_mergecols( GDanalysis.optcutoffpercs(TsI, curclass, curlabel), num2cell(repmat(mOptCutoffPercs, numel(TsI),1)));  
                        end 
                    end
                end
            end
            
            if MULTI.flag

                % MULTI-CLASS PERFORMANCE
                % =========================================================
                % Prepare for the OOT multi-class prediction by
                % concatenating prediction values for each sample 
                % across CV2 partitions
                for curlabel=1:MULTILABEL.dim
                    GDanalysis.grid.MultiSelNodeFreq(:,ll,curlabel) = GD.MultiGroupGridSelection{curlabel}.SelNodes;
                   
                    if (isfield(GRD,'NodeSelect') && ( GRD.NodeSelect.mode == 2 || GRD.NodeSelect.mode == 3 )) || combcell
                        MultiPred = []; 
                        MultiCV2Pred = []; 
                        MultiProb = []; 
                        MultiCV2Prob = cell(1,ngroups); 
                        if ~MULTI.BinBind
                            % Compute multi-group performance measures            
                            GDanalysis.multi_bestTR(f,d,curlabel) = nm_nanmean(GD.MultiGroupGridSelection{curlabel}.bestacc);
                            GDanalysis.multi_bestTS(f,d,curlabel) = nm_nanmean(GD.MultiGroupGridSelection{curlabel}.besttestparam);
                            % Store multi-group grid position
                            GDanalysis.multi_bestPpos{ll,curlabel} = GD.MultiGroupGridSelection{curlabel}.Npos;
                            % Select from multi-group prediction grid   
                            for zu=1:GD.MultiGroupGridSelection{curlabel}.Nodes
                                MultiPred = [MultiPred GD.MultiGroupGridSelection{curlabel}.bestpred{zu}];
                                MultiCV2Pred = [MultiCV2Pred GD.MultiGroupGridSelection{curlabel}.bestCV2pred{zu} ];
                                MultiProb = [ MultiProb GD.MultiGroupGridSelection{curlabel}.bestprob{zu} ];
                                for curclass=1:ngroups
                                    try
                                        MultiCV2Prob{curclass} = [ MultiCV2Prob{curclass} GD.MultiGroupGridSelection{curlabel}.bestCV2prob{zu,curclass}];
                                    catch
                                        fprintf('prob')
                                    end
                                end
                            end
                        else
                            GDanalysis.multi_bestTR(f,d,curlabel) = MultiBinBind.Mean_CVPerf;
                            GDanalysis.multi_bestTS(f,d,curlabel) = GD.MultiGroupGridSelection{curlabel}.besttestparam;
                            MultiPred = GD.MultiGroupGridSelection{curlabel}.bestpred;
                            MultiCV2Pred = GD.MultiGroupGridSelection{curlabel}.bestCV2pred;
                            MultiProb = GD.MultiGroupGridSelection{curlabel}.bestprob;
                            MultiCV2Prob = GD.MultiGroupGridSelection{curlabel}.bestCV2prob;
                        end
                        
                    else
                        %MultiProb = zeros([size(GD.MultiGroupGridSelection{curlabel}.bestprob{1}) GD.MultiGroupGridSelection{curlabel}.Nodes]);                        
                        if ~MULTI.BinBind
                            % Compute multi-group performance measures            
                            GDanalysis.multi_bestTR(f,d,curlabel) = GD.MultiGroupGridSelection{curlabel}.bestacc(1);
                            GDanalysis.multi_bestTS(f,d,curlabel) = GD.MultiGroupGridSelection{curlabel}.besttestparam(1);

                            % Select from multi-group prediction grid
                            MultiPred = GD.MultiGroupGridSelection{curlabel}.bestpred{1};
                            MultiCV2Pred = GD.MultiGroupGridSelection{curlabel}.bestCV2pred{1};
                            MultiProb = GD.MultiGroupGridSelection{curlabel}.bestprob{1};
                            MultiCV2Prob = GD.MultiGroupGridSelection{curlabel}.bestCV2prob;
                            % Store multi-group grid position
                            GDanalysis.multi_bestPpos(ll,curlabel) = GD.MultiGroupGridSelection{curlabel}.Npos(1);
                        else
                            GDanalysis.multi_bestTR(f,d,curlabel) = MultiBinBind.Mean_CVPerf;
                            GDanalysis.multi_bestTS(f,d,curlabel) = GD.MultiGroupGridSelection{curlabel}.besttestparam;
                            % Select multi-group prediction from binary optima
                            MultiPred = GD.MultiGroupGridSelection{curlabel}.bestpred;
                            MultiCV2Pred = GD.MultiGroupGridSelection{curlabel}.bestCV2pred;
                            MultiProb = GD.MultiGroupGridSelection{curlabel}.bestprob;
                            MultiCV2Prob = GD.MultiGroupGridSelection{curlabel}.bestCV2prob;
                        end
                    end
                    if ~RFE.CV2Class.EnsembleStrategy.AggregationLevel
                        MEnsDat = MultiPred; 
                    else
                        if iscell(MultiCV2Pred), MEnsDat = nk_cellcat(MultiCV2Pred,[],2); else, MEnsDat = MultiCV2Pred; end
                    end
                    % Concatenate multi-group prediction across CV2 perms
                    GDanalysis.multi_predictions(TsInd,curlabel) = cellmat_mergecols(GDanalysis.multi_predictions(TsInd,curlabel), num2cell(MEnsDat,2));
                    for g=1:ngroups
                        if ~RFE.CV2Class.EnsembleStrategy.AggregationLevel
                            MEnsDat = MultiProb(:,g,:); 
                        else
                            if iscell(MultiCV2Prob)
                                MEnsDat = MultiCV2Prob{g};
                            else
                                MEnsDat = MultiCV2Prob(:,g);
                            end
                        end
                        GDanalysis.multi_probabilities(TsInd, g, curlabel) = cellmat_mergecols(GDanalysis.multi_probabilities(TsInd, g, curlabel), num2cell(MEnsDat,2));
                    end
                    GDanalysis.grid.MultiCVPerf(:,ll,curlabel) = GD.MultiTR(:,curlabel);
                    GDanalysis.grid.MultiTSPerf(:,ll,curlabel) = GD.MultiTS(:,curlabel);
                    GDanalysis.grid.MultiERR_CVTSPerf(:,ll,curlabel) = GD.MultiERR(:,curlabel);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%% ENSEMBLE LEARNING %%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(ENSEMBLE) && (GDfl || ~batchflag)
            Meta = nk_GridSearchEnsemble(mapY, GD, k, g);
            oCVpathMeta = fullfile(tdir,[strout cvstr '_META_ID' id '.mat']);
            GDanalysis.METApaths{ll} = oCVpathMETA;
            fprintf('\nSaving %s', oCVpathMeta)
            if OCTAVE
                save(oCVpathMeta,'Meta', 'operm', 'ofold');
            else
                save(oCVpathMeta,'-v7.3','Meta', 'operm', 'ofold');
            end
            GDanalysis.Meta{ll} = Meta;
        end
        ll=ll+1;    
    end
end

if GDfl || batchflag == 0 || batchflag == 2

    GDanalysis.NumCV2Part = ol;

    % ********************* CV2 PERFORMANCE OVER GRID *********************
    GDanalysis.grid.se_CVPerf           = nm_nanstd(GDanalysis.grid.mean_CVPerf,3);
    GDanalysis.grid.se_TSPerf           = nm_nanstd(GDanalysis.grid.mean_TSPerf,3);
    GDanalysis.grid.se_Err_CVTSPerf     = nm_nanstd(GDanalysis.grid.mean_Err_CVTSPerf,3);
    GDanalysis.grid.se_Complexity       = nm_nanstd(GDanalysis.grid.mean_Complexity,3);
    GDanalysis.grid.se_CVDiversity      = nm_nanstd(GDanalysis.grid.mean_CVDiversity,3);
    GDanalysis.grid.se_TsDiversity      = nm_nanstd(GDanalysis.grid.mean_TsDiversity,3);
    GDanalysis.grid.mean_CVPerf         = nm_nanmean(GDanalysis.grid.mean_CVPerf,3);
    GDanalysis.grid.mean_TSPerf         = nm_nanmean(GDanalysis.grid.mean_TSPerf,3);
    GDanalysis.grid.mean_Err_CVTSPerf   = nm_nanmean(GDanalysis.grid.mean_Err_CVTSPerf,3);
    GDanalysis.grid.mean_Complexity     = nm_nanmean(GDanalysis.grid.mean_Complexity,3);
    GDanalysis.grid.mean_CVDiversity    = nm_nanmean(GDanalysis.grid.mean_CVDiversity,3);
    GDanalysis.grid.mean_TsDiversity    = nm_nanmean(GDanalysis.grid.mean_TsDiversity,3);
    GDanalysis.grid.SelNodeFreq         = bsxfun(@rdivide, nm_nansum(GDanalysis.grid.SelNodeFreq,3)*100, sum(nm_nansum(GDanalysis.grid.SelNodeFreq,3)));
    % Specifically account for the sequence optimizer                                     
    switch SVM.prog
        case 'SEQOPT'
            GDanalysis.grid.mean_SeqExamFreq    = nm_nanmean(GDanalysis.grid.mean_mSEQE,4);
            GDanalysis.grid.se_SeqExamFreq      = nm_nanmean(GDanalysis.grid.sd_mSEQE,4);
            GDanalysis.grid.mean_SeqGain        = nm_nanmean(GDanalysis.grid.mean_mSEQI,4);
            GDanalysis.grid.se_SeqGain          = nm_nanmean(GDanalysis.grid.sd_mSEQI,4);
            GDanalysis.grid.mean_SeqAbsUpper    = nm_nanmean(GDanalysis.grid.mean_mSEQAU,4);
            GDanalysis.grid.se_SeqAbsUpper      = nm_nanmean(GDanalysis.grid.sd_mSEQAU,4);
            GDanalysis.grid.mean_SeqAbsLower    = nm_nanmean(GDanalysis.grid.mean_mSEQAL,4);
            GDanalysis.grid.se_SeqAbsLower      = nm_nanmean(GDanalysis.grid.sd_mSEQAL,4);
            GDanalysis.grid.mean_SeqPercUpper   = nm_nanmean(GDanalysis.grid.mean_mSEQPU,4);
            GDanalysis.grid.se_SeqPercUpper     = nm_nanmean(GDanalysis.grid.sd_mSEQPU,4);
            GDanalysis.grid.mean_SeqPercLower   = nm_nanmean(GDanalysis.grid.mean_mSEQPL,4);
            GDanalysis.grid.se_SeqPercLower     = nm_nanmean(GDanalysis.grid.sd_mSEQPL,4);
            GDanalysis.grid.mean_SeqPerfGains   = nm_nanmean(GDanalysis.grid.mean_SeqPerfGains,4);
            [GDanalysis.CV2grid.caseprop_freq, GDanalysis.CV2grid.caseprop_node] = ...
                nk_ComputeEnsembleCasePropagationProbability(GDanalysis.caseprops,size(SVM.SEQOPT.C,2));
            DecValTrajConCat = cell2matpadnan(GDanalysis.decvaltraj);
            GDanalysis.CV2grid.decvaltraj       = nm_nanmedian(DecValTrajConCat, ndims(DecValTrajConCat) );
        case 'WBLCOX'
            GDanalysis.grid.mean_CutOffProb     = nm_nanmean(GDanalysis.grid.mean_mCutOffProb,3);
            GDanalysis.grid.se_CutOffProb       = nm_nanmean(GDanalysis.grid.mean_sdCutOffProb,3);
            GDanalysis.grid.mean_CutOffPerc     = nm_nanmean(GDanalysis.grid.mean_mCutOffPerc,3);
            GDanalysis.grid.se_CutOffPerc       = nm_nanmean(GDanalysis.grid.mean_sdCutOffPerc,3);
    end                                
                                        
    % This has to be changed to work in multi-label mode
    for h=1:nclass
        GDanalysis.best_CVperf{h} = mean(GDanalysis.bestTR{h}(GridUsed));
        GDanalysis.best_TSperf{h} = mean(GDanalysis.bestTS{h}(GridUsed));
        GDanalysis.best_Complexity{h} = mean(GDanalysis.bestComplexity{h}(GridUsed));
        GDanalysis.best_Error{h} = mean(GDanalysis.bestError{h}(GridUsed));
    end
    
    % ********************** ANALYSIS ACROSS PERMS ************************
    switch MODEFL
        case 'regression'
            lb = inp.labels;
            if MULTILABEL.flag && isfield(MULTILABEL,'sel')
               lb = inp.label(:,MULTILABEL.sel);
            end
            GDanalysis.Regr = nk_ComputeEnsembleProbability(GDanalysis.predictions(:,1,:), lb);

        case 'classification'
            
            labelh = zeros(size(label,1),nclass);
            %predh = zeros(size(label,1),nclass);
            for h=1:nclass
                % Build binary label vector
                if isfield(inp,'simFlag') && inp.simFlag
                    if numel(xCV.class{1,1}{h}.groups) == 2
                        ind1 = label == xCV.class{1,1}{h}.groups(1); ind2 = label == xCV.class{1,1}{h}.groups(2);
                        labelh(ind1,h) = 1; labelh(ind2,h) = -1;
                    else
                        ind1 = label == xCV.class{1,1}{h}.groups(1);
                        labelh(ind1,h) = 1; labelh(~ind1,h) = -1;
                    end
                else
                    if numel(CV.class{1,1}{h}.groups) == 2
                        ind1 = label == CV.class{1,1}{h}.groups(1); ind2 = label == CV.class{1,1}{h}.groups(2);
                        labelh(ind1,h) = 1; labelh(ind2,h) = -1;
                    else
                        ind1 = label == CV.class{1,1}{h}.groups(1);
                        labelh(ind1,h) = 1; labelh(~ind1,h) = -1;
                    end
                end

                labelhx = labelh(:,h); labelhx(labelh(:,h)<0)=0; Ix = find(labelh(:,h)); nIx = numel(Ix); 
                if ix>1
                    Px = GDanalysis.CV2grid.predictions; [~, Nx, ~] = size(Px); 
                    GDanalysis.CV2grid.mean_predictions(:,h) = nm_nanmean(Px(:,:,h),2);
                    GDanalysis.CV2grid.std_predictions(:,h)  = nm_nanstd(Px(:,:,h),2);
                    % I love anonymous functions - Compute performance measures
                    % for each permutation in the Grid.
                    try
                        %[ GDanalysis.CV2grid.Xsvm(:,:,h), GDanalysis.CV2grid.Ysvm(:,:,h) ] = arrayfun( @(j) perfcurve2(labelh(:,h), Px(:,j,h), 1), 1:Nx );
                        GDanalysis.CV2grid.CI_predictions(Ix,:,h)= cell2mat(arrayfun(@(i) percentile(Px(Ix(i),:,h),[2.5 97.5]),1:numel(Ix),'UniformOutput',false)');
                        GDanalysis.CV2grid.BAC(:,h)              = arrayfun( @(j) BAC( labelh(Ix,h),Px(Ix,j,h)-0.5), 1:Nx );
                        GDanalysis.CV2grid.sens(:,h)             = arrayfun( @(j) SENSITIVITY( labelh(Ix,h), Px(Ix,j,h)-0.5), 1:Nx );
                        GDanalysis.CV2grid.spec(:,h)             = arrayfun( @(j) SPECIFICITY( labelh(Ix,h), Px(Ix,j,h)-0.5), 1:Nx );
                        GDanalysis.CV2grid.accuracy(:,h)         = arrayfun( @(j) ACCURACY( labelh(Ix,h), Px(Ix,j,h)-0.5), 1:Nx );
                        GDanalysis.CV2grid.PPV(:,h)              = arrayfun( @(j) PPV( labelh(Ix,h), Px(Ix,j,h)-0.5 ), 1:Nx );
                        GDanalysis.CV2grid.NPV(:,h)              = arrayfun( @(j) NPV( labelh(Ix,h), Px(Ix,j,h)-0.5 ), 1:Nx );
                        GDanalysis.CV2grid.AUC(:,h)              = arrayfun( @(j) fastAUC( labelhx(Ix), Px(Ix,j,h)-0.5, 1), 1:Nx );
                        GDanalysis.CV2grid.DOR(:,h)              = arrayfun( @(j) DOR( labelh(Ix,h), Px(Ix,j,h)-0.5 ), 1:Nx );
                    catch
                        warning('CVdatamats of more than one CV2 permutation are needed.')
                    end
                end
                switch SVM.prog
                    case 'WBLCOX'
                        GDanalysis.BinClass{h} = nk_ComputeEnsembleProbability(GDanalysis.predictions(:,h), labelh(:,h),[], GDanalysis.CV2cutoffs, GDanalysis.optcutoffpercs);
                        GDanalysis.BinClass{h}.Time.mean_predicted_times = cell2matpadnan(GDanalysis.predtimes(Ix,h));
                        if size(GDanalysis.BinClass{h}.Time.mean_predicted_times,2)==1, GDanalysis.BinClass{h}.Time.mean_predicted_times=[GDanalysis.BinClass{h}.Time.mean_predicted_times nan(nIx,1)]; end
                        GDanalysis.BinClass{h}.Time.mean_predicted_times = nm_nanmedian(GDanalysis.BinClass{h}.Time.mean_predicted_times,2); 
                        GDanalysis.BinClass{h}.Time.observed_times = time2event(Ix);
                        ind = GDanalysis.BinClass{h}.Time.mean_predicted_times>0;
                        [rr, pp, rl, ru]                       = corrcoef(GDanalysis.BinClass{h}.Time.observed_times(ind), GDanalysis.BinClass{h}.Time.mean_predicted_times(ind) );
                        GDanalysis.BinClass{h}.Time.R2         = SCC(GDanalysis.BinClass{h}.Time.observed_times(ind), GDanalysis.BinClass{h}.Time.mean_predicted_times(ind));
                        GDanalysis.BinClass{h}.Time.r          = rr(1,2);
                        GDanalysis.BinClass{h}.Time.p          = pp(1,2);
                        GDanalysis.BinClass{h}.Time.r_95CI_low = rl(1,2);
                        GDanalysis.BinClass{h}.Time.r_95CI_up  = ru(1,2);
                        %% Convert correlation coefficient to T value
                        GDanalysis.BinClass{h}.Time.t          = rr(1,2) * sqrt((sum(ind)-2) / (1 - rr(1,2)*rr(1,2)));
                        
                    otherwise
                        % Do performance stats over the entire experiment
                        GDanalysis.BinClass{h} = nk_ComputeEnsembleProbability(GDanalysis.predictions(:,h), labelh(:,h));
                end
            end
    end

    if MULTI.flag && nclass > 1

        % ****************** MAX SELECTION OVER CV2 GRID ******************
        TRvec = GDanalysis.multi_bestTR(GridUsed);
        TSvec = GDanalysis.multi_bestTS(GridUsed);
        GDanalysis.best_MultiCVperf = nm_nanmean(TRvec(:));
        GDanalysis.best_MultiTSperf = nm_nanmean(TSvec(:));
        GDanalysis.best_sdMultiCVperf = nm_nanstd(TRvec(:));
        GDanalysis.best_sdMultiTSperf = nm_nanstd(TSvec(:));

        % ************** MULTI-CLASS PERF across CV2-PERMS ****************
        %multi_pred = GDanalysis.multi_predictions;
        multi_CV2pred = GDanalysis.multi_predictions;
        multi_CV2prob = GDanalysis.multi_probabilities;
        % ** Mean multi-classification performance across CV2 partitions **
        GDanalysis.grid.seMultiCVPerf       = nm_nanstd(GDanalysis.grid.MultiCVPerf,2);
        GDanalysis.grid.seMultiTSPerf       = nm_nanstd(GDanalysis.grid.MultiTSPerf,2);
        GDanalysis.grid.seMultiERR_CVTSPerf = nm_nanstd(GDanalysis.grid.MultiERR_CVTSPerf,2);
        GDanalysis.grid.seMultiCVDiversity  = nm_nanstd(GDanalysis.grid.MultiCVDiversity,2);
        GDanalysis.grid.seMultiTsDiversity  = nm_nanstd(GDanalysis.grid.MultiTsDiversity,2);
        GDanalysis.grid.seMultiComplexity   = nm_nanstd(GDanalysis.grid.mean_Complexity,2);
        GDanalysis.grid.MultiCVPerf         = nm_nanmean(GDanalysis.grid.MultiCVPerf,2);
        GDanalysis.grid.MultiTSPerf         = nm_nanmean(GDanalysis.grid.MultiTSPerf,2);
        GDanalysis.grid.MultiERR_CVTSPerf   = nm_nanmean(GDanalysis.grid.MultiERR_CVTSPerf,2);
        GDanalysis.grid.MultiCVDiversity    = nm_nanmean(GDanalysis.grid.MultiCVDiversity,2);
        GDanalysis.grid.MultiTsDiversity    = nm_nanmean(GDanalysis.grid.MultiTsDiversity,2);
        GDanalysis.grid.MultiComplexity     = nm_nanmean(GDanalysis.grid.mean_Complexity,2);
        GDanalysis.grid.MultiSelNodeFreq    = nm_nanmean(GDanalysis.grid.MultiSelNodeFreq,2);
        % Convert OOT predictions to probabilities for class membership
        GDanalysis = nk_MultiPerfComp(GDanalysis, multi_CV2pred, label, ngroups);
        tProbCat = nk_cellcat(multi_CV2prob(:,1),[],1);
        mProb = zeros([lx size(tProbCat,2) ngroups]);
        Ix = ~cellfun(@isempty, multi_CV2prob(:,1));
        for g=1:ngroups
            mProb(Ix,:,g) =  nk_cellcat(multi_CV2prob(:,g),[],1);
        end
        mProbI = squeeze(nm_nanmean(mProb,2)); mProbI(~Ix,:)=NaN;
        GDanalysis = nk_MultiPerfComp(GDanalysis, mProbI, label, ngroups, 'prob');
    end

    % **************************** SAVE DATA *****************************
    if ~batchflag
        fprintf('\nSaving %s', GDanalpth);
        GDanalysis.path     = GDanalname;
        if OCTAVE
            save(GDanalpth,'GDanalysis');
        else
            save(GDanalpth,'-v7.3','GDanalysis');
        end
    end
else
    if fromData
        xNM.simanalysis = GDanalysis;
    else 
        GDanalysis = [];
    end
end
