% =========================================================================
% FORMAT [InputParam, TrainedParam] = nk_PerfPreprocess_core(InputParam, 
%                                                             TrParam, act)
% =========================================================================
% This is the main preprocessing module of NeuroMiner. It orchestrates the
% sequence of pre-processing steps to be applied to the data, collects and 
% organizes both the preprocessed data and the computed preprocessing 
% parameters, and manages CV and OOCV modes of operation. If hyperparameters 
% are part of the preprocessing, the function manages the creation of data 
% shelves where different version of data processed using the given 
% hyperparameter combination are stored. This results in a tree structure 
% of preprocessed data and respective parameters.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2024

function [SrcParam, InputParam, TrParam] = ...
              nk_PerfPreprocessObj_core(SrcParam, InputParam, TrParam, act)

global VERBOSE TEMPL SVM MODEFL SYNTH 

%% Prepare and check input params
paramfl = false; tsfl = false; trfl = false; cfl = false;

% Check if training params are empty and therefore have to be computed.
if ~exist('TrParam','var') || isempty(TrParam)
    if VERBOSE; fprintf('\nNo Training parameters found. Training mode enabled.'); end
    TrParam = []; 
else
    fprintf('\nTraining parameters found. Out-of-sample mode enabled.');   
    paramfl = true; 
end

% Is therey any data to be processed?
if ~exist('InputParam','var') || isempty(InputParam), return; end

% Is training data in the input parameters?
if isfield(InputParam,'Tr'), trfl=true; end

% Are multiple test data sets defined?
nTs = 0; if isfield(InputParam,'Ts'), tsfl=true; if iscell(InputParam.Ts), nTs = size(InputParam.Ts,2); end; end

% Has calibration data been transferred? [not used currently]
nC = 0; 
if isfield(InputParam,'C')
   
    cfl=true; 
    if iscell(InputParam.C)
        nC = size(InputParam.C,2); 
    end
end

if ~trfl && ~tsfl, error('\nNo Preprocessing performed because training and test data are missing.\nCheck your parameters!!!'); end

% Has an action sequence been provided? 
if iscell(act), nact = length(act); else, nact = 1; acti = act; end

% Do the training parameters and does the training parameters sequence
% match action sequence?
if iscell(TrParam), nTrParam = length(TrParam); else, nTrParam = 1; TrParami = []; end
if nact ~= nTrParam && paramfl, if VERBOSE; fprintf('\nLength of action sequence has to match training parameters. Abort!!!'); end; return; end

% Do we need to deal with multiple labels?
switch MODEFL
    case 'classification'
        nTrL = size(SrcParam.MultiTrainLabel,2);
    case 'regression'
        nTrL = size(SrcParam.TrainLabel,2);
end
% if we deal with multiple labels, do we need to create one processed data version
% for each label or can we use a single processed data version for all
% labels, which saves time and memory. For this, we need to find out if
% there are label interactions during preprocessing. Label interactions are
% defined in the parent function nk_GenPreprocSequence.m
if nTrL> 1 && nk_Check4LabelInteraction(InputParam.P) 
    nL = nTrL; 
    if VERBOSE, fprintf('\nMulti-label mode detected: Processing %g labels', nL); end
else
    nL = 1;
end

% Do we have to manage synthetic data?
synthfl= false; 
if ~isempty(SYNTH) && SYNTH.flag == 1 && any(SrcParam.synth_activated) 
    synthfl= true; 
elseif isfield(SVM,'ADASYN') && isfield(SVM.ADASYN, 'flag') && SVM.ADASYN.flag == 1 && any(SrcParam.synth_activated) 
    synthfl= true; 
end

% done with checking and preparing, now let's do the job
% -------------------------------------------------------------------------

if nact>1, fprintf('\t...Execute preprocessing sequence: '); end

if isfield(InputParam, 'Ts') && ~isempty(InputParam.Ts) && isequaln(InputParam.Tr, InputParam.Ts{1})
    copy_ts = true;
else
    copy_ts = false;
end

ActParam = struct('trfl', trfl, ...
                  'paramfl', paramfl, ...
                  'tsfl', tsfl, ...
                  'cfl', cfl, ...
                  'nTs', nTs', ...
                  'nC', nC, ...
                  'i',0, ...
                  'synthfl', synthfl, ...
                  'curlabel', 1, ...
                  'label_interaction', nL>1, ...
                  'copy_ts', copy_ts);

% if writeout is set to yes the function will write out data and parameters
% at each hyperparameter combination. This is a debugging mode for the
% NM developers.
writeout = false;

for i=1:nact
    
    ActParam.i=i;
    tStart = tic; 
    
    if iscell(act), acti = act{i}; end
    % Out-of-sample mode: get current training parameters
    if iscell(TrParam) && paramfl, TrParami = TrParam{i}; end
    
    tic; if VERBOSE;fprintf('\nStep #%g/%g:',i,nact); end
    funcstr = sprintf('act_%s', acti); 
    
    % Check whether multiple shelfs of data are available
    if iscell(InputParam.Tr) 
        
        % Check for parameter range at current step
        nP = size(InputParam.Tr,1); O = []; nO = 1;
        if isfield(InputParam.P{i},'opt') && ~isempty(InputParam.P{i}.opt)
            O = InputParam.P{i}.opt;
            nO = size(InputParam.P{i}.opt,1);
            if VERBOSE && nO>1
                fprintf('%g hyperparameters detected at preprocessing step %g.',nO, i); 
            end
        end
        
        % Combinations of shelves and current parameter range
        nComb = nP * nO; 
        Tr = cell(nComb,nL); 
        if isfield(InputParam,'Yw'), Yw = cell(nComb,nL); end
        if isfield(InputParam,'C'), C = cell(nComb,nL); end
        if nTs>0, Ts = cell(nComb,nTs,nL); end
        TrParamij = cell(nComb,nL); ll = 1;
        
        % Prepare Input parameter container for multi-parameter operations
        InputParamj.P               = InputParam.P;
        InputParamj.CV1perm         = InputParam.CV1perm;
        InputParamj.CV1fold         = InputParam.CV1fold;
        InputParamj.curclass        = InputParam.curclass;
        if nTs > 1, InputParamj.Ts  = cell(1,nTs);   end
        
        % Do we have to take synthetic data generation into account?
        if synthfl
            TrainLabelSyn = cell(nComb,nL); 
            CovarSyn = cell(nComb,nL);
        end
        
        % Loop through multiple labels if needed
        for nl=1:nL
            
            kl=ll;
            ActParam.curlabel = nl;
            
            % Loop through parameter space
            for j=1:nP
                
                ActParam.j = j;
                % Get training and test data from current parameter shelf
                jj=1; if size(InputParam.Tr,1) > 1; jj = j; end
                if synthfl&& i == 1
                   InputParamj.Tr = [ InputParam.Tr{jj}; InputParam.TrSyn{jj} ]; 
                else 
                   if i==1    
                        InputParamj.Tr = InputParam.Tr{jj}; 
                   else
                        InputParamj.Tr = InputParam.Tr{jj,nl}; 
                   end
                end
                if isfield(InputParam,'Yw') 
                    try
                        if i==1
                            InputParamj.Yw = InputParam.Yw{j}; 
                        else
                            InputParamj.Yw = InputParam.Yw{j,nl}; 
                        end
                    catch
                        if i==1
                            InputParamj.Yw = InputParam.Yw;
                        else
                            InputParamj.Yw = InputParam.Yw{nl};
                        end
                    end
                end
                if ~isempty(TEMPL) && isfield(TEMPL,'Param') && isfield(TEMPL,'execute_template_creation') && ~TEMPL.execute_template_creation
                    if isstruct(TEMPL.Param(InputParam.curclass).TrainedParam{i})
                        ActParam.Templ = TEMPL.Param(InputParam.curclass).TrainedParam{i};
                    else
                        ActParam.Templ = TEMPL.Param(InputParam.curclass).TrainedParam{i}{j}; 
                    end
                end
                if nTs>1
                    for k = 1:nTs
                        if size( InputParam.Ts,1)==1
                            if i==1
                                InputParamj.Ts{k} = InputParam.Ts{1,k};
                            else
                                InputParamj.Ts{k} = InputParam.Ts{1,k,nl};
                            end
                        else
                            if i==1
                                InputParamj.Ts{k} = InputParam.Ts{j,k};
                            else
                                InputParamj.Ts{k} = InputParam.Ts{j,k, nl};
                            end
                        end
                    end
                elseif nTs==1
                    if size( InputParam.Ts,1 )==1
                        if i==1
                            InputParamj.Ts = InputParam.Ts{1};
                        else
                            InputParamj.Ts = InputParam.Ts{1,nl};
                        end
                    else
                        if i==1
                            InputParamj.Ts = InputParam.Ts{j};
                        else
                            InputParamj.Ts = InputParam.Ts{j, nl};
                        end
                    end
                end
                if isfield(InputParam,'C')
                    InputParam.C = InputParam.C;
                    if i==1    
                        InputParamj.C = InputParam.C{jj}; 
                    else
                        InputParamj.C = InputParam.C{jj,nl}; 
                   end
                    
                end
                % ... and now apply multiple parameters to shelf thus creating
                % multi-shelf versions of the data 
                for l = 1: nO
                    if ~isempty(O)
                        ActParam.opt = O(l,:);
                    elseif isfield(ActParam,'opt') 
                        ActParam = rmfield(ActParam,'opt');
                    end
                    if paramfl, llTrParam = TrParami{kl,nl}; else, llTrParam = []; end
                    [ SrcParam, Out, TrParamij{kl,nl}, ActParam ] = feval( funcstr, SrcParam, InputParamj, TrParam, llTrParam, ActParam );
                    
                    Tr{kl,nl} = Out.Tr;
                    if writeout
                        filepath = fullfile(pwd, sprintf('PreprocDataMat_CV2-%g-%g_CV1-%g-%g_Step%g_nP%g_HypP%g.mat', ...
                            CVPOS.CV2p, CVPOS.CV2f, CVPOS.CV1p, CVPOS.CV1f, i, j, l));
                        fprintf('\nSaving preprocessing data:\n%s', filepath);
                        save(filepath, "Out")
                    end
                    if synthfl, TrainLabelSyn{kl,nl} = SrcParam.TrainLabelSyn{ActParam.j}; CovarSyn{kl,nl} = SrcParam.covarsSyn{ActParam.j}; end
                    if isfield(Out,'Yw'), Yw{kl,nl} = Out.Yw; end
                    if isfield(Out, 'C'), C{kl,nl} = Out.C; end
                    if nTs>1
                        for k = 1:nTs
                            Ts{kl,k,nl} = Out.Ts{k};
                        end
                    elseif nTs==1
                        Ts{kl,nl} = Out.Ts;
                    end
                    kl=kl+1;
                end
            end
        end
        InputParam.Tr = Tr; 
        if synthfl, SrcParam.TrainLabelSyn = TrainLabelSyn; SrcParam.covarsSyn = CovarSyn; end
        if isfield(InputParam,'Yw'), InputParam.Yw = Yw; end
        if isfield(InputParam,'C'), InputParam.C = C; end

        if nTs > 0; InputParam.Ts = Ts; end
        TrParami = TrParamij;
    
    % DO WE HAVE A PARAMETER RANGE FOR THE CURRENT STEP?    
    elseif isfield(InputParam.P{i},'opt') && ~isempty(InputParam.P{i}.opt)
        O = InputParam.P{i}.opt; % this is the parameter range to be processed
        nO = size(InputParam.P{i}.opt,1);
        Tr = cell(nO,nL); 
        % Is there are weighting vector?
        if isfield(InputParam,'Yw'), Yw = cell(nO,nL); end
        if isfield(InputParam,'C'), C = cell(nO,nL); end
        if nTs > 0, Ts = cell(nO,nTs,nL); end

        TrParamij = cell(nO,nL); ll = 1;
        %Do we have to use Adasyn ?
        if synthfl
            % Prepare containers for synthetic label and covars expansion
            TrainLabelSyn = cell(nO,nL); 
            CovarSyn = cell(nO,nL);
            if i==1
                % Concatenate original and synthetic training data
                InputParam.Tr = [InputParam.Tr; InputParam.TrSyn];
            end
        end
        % -----------------------------------------------------------------
        % Loop through labels if necessary
        for nl = 1:nL
            ActParam.curlabel = nl;
            kl=ll;
            % Now, loop through parameter range at given preprocessing step
            for l = 1: nO
                % PREPARE PARAMS
                ActParam.opt = O(l,:); % Retrieve current processing params
                % Is there any template parameter structure?
                if ~isempty(TEMPL) && isfield(TEMPL,'Param') && isfield(TEMPL,'execute_template_creation') && ~TEMPL.execute_template_creation
                    ActParam.Templ = TEMPL.Param(InputParam.curclass).TrainedParam{i}; 
                end
                % Are there any precomputed parameters?
                if paramfl, llTrParam = TrParami{kl,nl}; else, llTrParam = []; end
                % -----------------------------------------------------------------------------------------------------------------
                % PROCESS: Pass everything to processing step (funcstr)
                [ SrcParam, Out, TrParamij{kl,nl}, ActParam ] = feval( funcstr, SrcParam, InputParam, TrParam, llTrParam, ActParam );
                % -----------------------------------------------------------------------------------------------------------------
                if writeout
                    filepath = fullfile(pwd,sprintf('PreprocDataMat_CV2-%g-%g_CV1-%g-%g_Step%g_HypP%g.mat', ...
                        CVPOS.CV2p, CVPOS.CV2f, CVPOS.CV1p, CVPOS.CV1f, i, l));
                    fprintf('\nSaving preprocessing data:\n%s', filepath);
                    save(filepath, "Out")
                end
                % POST-PROCESS
                TrParamij{kl,nl}.i_opt = O(l,:); % Save in parameter shelf
                Tr{kl,nl} = Out.Tr; % Save in data shelf
                % Take care of adasyn labels and covars expansion if needed.
                if synthfl
                    TrainLabelSyn{kl, nl} = SrcParam.TrainLabelSyn; 
                    if ~isempty(SrcParam.covars)
                        CovarSyn{kl, nl} = SrcParam.covarsSyn; 
                    end
                end
                % Take care of weighting vector to be also stored in shelf
                if isfield(Out,'Yw'), Yw{kl,nl} = Out.Yw; end
                if isfield(Out,'C'), C{kl,nl} = Out.C; end
                % Now transfer proper training, CV1 test and CV2 test data to
                % data shelves
                if nTs>1
                    for k = 1:nTs, Ts{kl,k,nl} = Out.Ts{k}; end
                elseif nTs==1
                    Ts{kl,1,nl} = Out.Ts;
                end
                kl=kl+1;
                % -------------------------------------------------------------
            end
        end
        % Save processed training data shelves back in InputParam.Tr, which is now a cell array 
        InputParam.Tr = Tr; 
        % Adasyn: Store expanded labels and covars back in SrcParam
        if synthfl, SrcParam.TrainLabelSyn = TrainLabelSyn; SrcParam.covarsSyn = CovarSyn; end
        % Save CV1 and CV2 test data in shelves
        if isfield(InputParam,'Yw'), InputParam.Yw = Yw; end
        if nC > 0, InputParam.C = C; end
        if nTs > 0; InputParam.Ts = Ts; end
        TrParami = TrParamij; % Save processing parameters of curent parameters range loop
    else
        if synthfl&& i==1
            InputParam.Tr = [InputParam.Tr; InputParam.TrSyn];
        end
        if ~isempty(TEMPL) && isfield(TEMPL,'Param') && isfield(TEMPL,'execute_template_creation') && ~TEMPL.execute_template_creation
            ActParam.Templ = TEMPL.Param(InputParam.curclass).TrainedParam{i}; 
        end
        [ SrcParam, InputParam, TrParami, ActParam ] = feval( funcstr, SrcParam, InputParam, TrParam, TrParami, ActParam );
        if writeout
            Out = InputParam;
            filepath = fullfile(pwd,sprintf('PreprocDataMat_CV2-%g-%g_CV1-%g-%g_Step%g.mat', ...
                CVPOS.CV2p, CVPOS.CV2f, CVPOS.CV1p, CVPOS.CV1f, i));
            fprintf('\nSaving preprocessing data:\n%s', filepath);
            save(filepath, "Out");
            clear Out;
        end
    end
    
    % No Out-of-sample mode => build up TrParam structure array
    if ~paramfl, TrParam{i} = TrParami; end; TrParami = []; 
    
    if VERBOSE, fprintf('\tDone in %g seconds.', toc(tStart)); else, fprintf('.'); end
    
end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam] = act_impute(SrcParam, InputParam, ~, TrParami, actparam)
%global VERBOSE
trfl     = actparam.trfl;
tsfl     = actparam.tsfl;
copy_ts  = actparam.copy_ts;
%paramfl  = actparam.paramfl;
i        = actparam.i;
tsproc   = false;
InputParam.P{i}.IMPUTE.X = InputParam.Tr;

if trfl 
    [InputParam.Tr, TrParami] = nk_PerfImputeObj(InputParam.Tr, InputParam.P{i}.IMPUTE); 
    TrParami.X = InputParam.Tr;
    if tsfl, tsproc = true; end  
end
if tsproc 
    if copy_ts
        Ts = nk_PerfImputeObj(InputParam.Ts(2:end), TrParami); 
        InputParam.Ts{1} = InputParam.Tr;
        InputParam.Ts(2:end) = Ts;
    else
        InputParam.Ts = nk_PerfImputeObj(InputParam.Ts, TrParami); 
    end
end
end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam] = act_labelimpute(SrcParam, InputParam, ~, ~, actparam)
global MODEFL
i  = actparam.i;
InputParam.P{i}.LABELIMPUTE.X = InputParam.Tr;
switch MODEFL
    case 'classification'
        L = SrcParam.MultiTrainLabel;
    case 'regression'
        L = SrcParam.TrainLabel;
end
if actparam.synthfl
    L = [L; SrcParam.TrainLabelSyn{actparam.j}];
end
[SrcParam.TrL_imputed, TrParami] = nk_PerfImputeLabelObj(L, InputParam.P{i}.LABELIMPUTE); 

end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam ] = act_scale(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;
if ~isfield(InputParam.P{i}.SCALE,'AcMatFl'), InputParam.P{i}.SCALE.AcMatFl = false; end
if ~isfield(InputParam.P{i}.SCALE,'ZeroOne'), InputParam.P{i}.SCALE.ZeroOne = 1;     end
if ~isfield(InputParam.P{i}.SCALE,'zerooutflag'), InputParam.P{i}.SCALE.zerooutflag = 1;     end
if VERBOSE
    if InputParam.P{i}.SCALE.ZeroOne == 1, rangestr = '[0, 1]'; else, rangestr = '[-1, 1]'; end
    if InputParam.P{i}.SCALE.AcMatFl, fprintf('\tScaling whole matrix to %s ...', rangestr); else, fprintf('\tScaling each feature to %s ...', rangestr); end
end
if paramfl && tsfl && isfield(TrParami,'minY') && ...
        isfield(TrParami,'maxY'), tsproc = true;
else
    if trfl, [InputParam.Tr, TrParami] = nk_PerfScaleObj(InputParam.Tr, InputParam.P{i}.SCALE); end
    if tsfl, tsproc = true; end  
end
if tsproc, InputParam.Ts = nk_PerfScaleObj(InputParam.Ts, TrParami); end
if isfield(InputParam, 'C') %CALIB.calibflag && CALIB.preprocstep > i
    InputParam.C = nk_PerfScaleObj(InputParam.C, InputParam.P{i}.SCALE);
end
end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam ] = act_normalize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc = false;

if VERBOSE; fprintf('\tNormalizing to group mean(s) ...'); end
            
if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(InputParam.P{i},'TsInd')
    tsproc = true;
elseif trfl
    if actparam.synthfl&& isfield(InputParam.P{i},'TrInd') && ~isempty(InputParam.P{i}.TrInd)
        if iscell(SrcParam.covarsSyn)
            InputParam.P{i}.TrInd = [InputParam.P{i}.TrInd; SrcParam.covarsSyn{actparam.j}(:,InputParam.P{i}.IND)]; 
        else
            InputParam.P{i}.TrInd = [InputParam.P{i}.TrInd; SrcParam.covarsSyn(:,InputParam.P{i}.IND)]; 
        end
    end
    [InputParam.Tr, TrParami] = nk_PerfNormObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end 
else
    if VERBOSE; fprintf('\n not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfNormObj(InputParam.Ts, TrParami); end

end
% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_copyfeat(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

if VERBOSE;fprintf('\tCopy features ...'); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_unitnormalize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE; fprintf('\tNormalizing to unit vector ...'); end 

if paramfl && tsfl && isfield(TrParami,'normY') 
    tsproc = true;
elseif isfield(InputParam.P{i},'METHOD')
    if trfl, [InputParam.Tr, TrParami] = nk_PerfUnitNormObj(InputParam.Tr, InputParam.P{i}); end
    if tsfl, tsproc = true; end 
else
    if VERBOSE;fprintf('\n not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfUnitNormObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_remmeandiff(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE CALIB
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE; fprintf('\tCorrect groups for offsets from global mean(s) ...'); end

if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(TrParami,'meanG')
    tsproc = true;
elseif trfl
     if actparam.synthfl&& isfield(InputParam.P{i},'sTrInd') && ~isempty(InputParam.P{i}.sTrInd)
        if iscell(SrcParam.covarsSyn)
            InputParam.P{i}.sTrInd = [InputParam.P{i}.sTrInd; SrcParam.covarsSyn{actparam.j}(:,InputParam.P{i}.sIND)]; 
        else
            InputParam.P{i}.sTrInd = [InputParam.P{i}.sTrInd; SrcParam.covarsSyn(:,InputParam.P{i}.sIND)]; 
        end
     end
    [InputParam.Tr, TrParami] = nk_PerfRemMeanDiffObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end 
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfRemMeanDiffObj(InputParam.Ts, TrParami); end
if isfield(InputParam,'C') && CALIB.calibflag % && CALIB.preprocstep > i
    INcalib = InputParam.P{i};
    INcalib.sTrInd = SrcParam.covars_cocv(:,INcalib.sIND);
    INcalib.dTrInd = SrcParam.covars_cocv(:,INcalib.dIND);
    InputParam.C = nk_PerfRemMeanDiffObj(InputParam.C, INcalib);
end
end

% =================================================================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_discretize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE;
    fprintf('\tDiscretizing to range [ %g : %g : %g ] ...', ...
        InputParam.P{i}.DISCRET.binstart, ...
        InputParam.P{i}.DISCRET.binsteps, ...
        InputParam.P{i}.DISCRET.binstop)
end
if paramfl && tsfl && isfield(TrParami,'mY') && isfield(TrParami,'sY') 
    tsproc = true;
elseif trfl
    [InputParam.Tr, TrParami] = nk_PerfDiscretizeObj(InputParam.Tr, InputParam.P{i});
    if tsfl, tsproc = true; end
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfDiscretizeObj(InputParam.Ts, TrParami); end

end

% =================================================================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_symbolize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
nTs     = actparam.nTs;
tsproc = false;    

if VERBOSE; fprintf('\tSymbolizing ...'); end

if paramfl && tsfl && isfield(TrParami,'sBin')
    tsproc=true;
elseif trfl
    [ InputParam.Tr, TrParami ] = nk_PerfDiscretizeObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfDiscretizeObj(InputParam.Ts, TrParami); end

end   

% =================================================================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_standardize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc = false;            

if isfield(actparam,'opt')
    InputParam.P{i}.WINSOPT = actparam.opt;
end

if VERBOSE; fprintf('\tStandardizing data ...'); end
if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(TrParami,'stdY')
    tsproc = true;
else
    if trfl
        % Check for ADASYN
        if actparam.synthfl&& isfield(InputParam.P{i},'sTrInd') && ~isempty(InputParam.P{i}.sTrInd)
            if iscell(SrcParam.covarsSyn)
                InputParam.P{i}.sTrInd = [InputParam.P{i}.sTrInd; SrcParam.covarsSyn{actparam.j}(:,InputParam.P{i}.sIND)]; 
            else
                InputParam.P{i}.sTrInd = [InputParam.P{i}.sTrInd; SrcParam.covarsSyn(:,InputParam.P{i}.sIND)]; 
            end
        end
        [InputParam.Tr, TrParami] = nk_PerfStandardizeObj(InputParam.Tr, InputParam.P{i});
    end
    if tsfl, tsproc = true; end
end
if tsproc, InputParam.Ts = nk_PerfStandardizeObj(InputParam.Ts, TrParami); end

end

% =================================================================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_reducedim(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
copy_ts  = actparam.copy_ts;
tsproc  = false;

if isfield(actparam,'DIMOPT') && ~isempty(actparam.DIMOPT)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DR.opt = actparam.DIMOPT;
elseif isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DR.opt = actparam.opt;
end 

% We need the template parameters for the optional Procrustes rotation
if isfield(actparam,'Templ') 
    if isfield(actparam.Templ,'mpp')
        InputParam.P{i}.DR.mpp_template = actparam.Templ.mpp;
    elseif iscell(actparam.Templ) && isfield(actparam.Templ{1},'mpp')
        InputParam.P{i}.DR.mpp_template = actparam.Templ{1}.mpp;
    end
end

if isfield(actparam,'PRUNE'), actparam = rmfield(actparam,'PRUNE'); end

InputParam.P{i}.TrX = InputParam.Tr;

if VERBOSE    
    if trfl
        featnum = size(InputParam.Tr,2); 
    elseif tsfl
        if iscell(InputParam.Ts)
            featnum = size(InputParam.Ts{1},2); 
        else
            featnum = size(InputParam.Ts,2); 
        end
    end
    fprintf('\tDimensionality reduction [%s]: %g input features -> ', ...
                                    InputParam.P{i}.DR.RedMode, featnum); 
end

DIMRED = InputParam.P{i};

if paramfl && tsfl && isfield(InputParam.P{i},'DR') && ...
        isfield(TrParami,'mpp')
    tsproc = true;
elseif trfl
     if actparam.synthfl
        if iscell(SrcParam.TrainLabelSyn)
            DIMRED.DR.labels = [InputParam.P{i}.DR.labels(:,actparam.curlabel); SrcParam.TrainLabelSyn{actparam.j}(:,actparam.curlabel)]; 
        else
            DIMRED.DR.labels = [InputParam.P{i}.DR.labels(:,actparam.curlabel); SrcParam.TrainLabelSyn(:,actparam.curlabel)]; 
        end
     else
         DIMRED.DR.labels = InputParam.P{i}.DR.labels(:,actparam.curlabel); 
     end
    [InputParam.Tr, TrParami] = nk_PerfRedObj(InputParam.Tr, DIMRED); 
    if tsfl,tsproc = true; end   
end 
TrParami.mappedX = InputParam.Tr;

if tsproc
    if copy_ts
        InputParam.Ts{1} = InputParam.Tr;
        Ts = nk_PerfRedObj(InputParam.Ts(2:end), TrParami); 
        InputParam.Ts(2:end) = Ts;
    else
        InputParam.Ts = nk_PerfRedObj(InputParam.Ts, TrParami); 
    end
end

end

% =================================================================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_extdim(SrcParam, InputParam, TrParam, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
copy_ts  = actparam.copy_ts;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.opt = actparam.opt;
end 

% Extract weight vector from last TrParam
if isfield(actparam,'j'); j = actparam.j; end
if exist('j','var')
    if isfield(TrParam{i-1}{j},'mpp')
        InputParam.P{i}.mpp = TrParam{i-1}{j}.mpp;
    end
else
    if isfield(TrParam{i-1},'mpp')
        InputParam.P{i}.mpp = TrParam{i-1}.mpp;
    end
end

if VERBOSE    
    if trfl
        featnum = size(InputParam.Tr,2); 
    elseif tsfl
        if iscell(InputParam.Ts)
            featnum = size(InputParam.Ts{1},2); 
        else
            featnum = size(InputParam.Ts,2); 
        end
    end
    fprintf('\tExtraction of subspaces from reduced data: %g input features -> ', featnum); 
end

if paramfl && tsfl && isfield(TrParami,'indNonRem')
    tsproc = true;
elseif trfl 
    [InputParam.Tr, TrParami] = nk_PerfExtDimObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl,tsproc = true; end   
end 

if tsproc 
    if copy_ts
        InputParam.Ts{1} = InputParam.Tr;
        Ts = nk_PerfExtDimObj(InputParam.Ts(2:end), TrParami); 
        InputParam.Ts(2:end)=Ts;
    else
        InputParam.Ts = nk_PerfExtDimObj(InputParam.Ts, TrParami); 
    end
end
end

% =================================================================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_correctnuis(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE CALIB 

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
copy_ts = actparam.copy_ts;
tsproc  = false;            

IN = InputParam.P{i}; 

switch IN.METHOD
    case 1
        % METHOD 1: PARTIAL CORRELATIONS.
        IN.revertflag = 0; if isfield(InputParam.P{i},'COVDIR') && InputParam.P{i}.COVDIR ==1, IN.revertflag = 1; end
        IN.nointercept = 1; if isfield(InputParam.P{i},'INTERCEPT'), IN.nointercept = ~InputParam.P{i}.INTERCEPT; end
        if VERBOSE
            fprintf('\tPartial correlations ...'); 
            if IN.revertflag, fprintf(' introduce covariate effects ...'); else, fprintf(' remove covariate effects ...'); end
            if IN.nointercept, fprintf(' intecept disabled ...'); else, fprintf(' intercept enabled ...'); end
        end
        if isfield(InputParam.P{i},'BETAEXT')
            if VERBOSE;fprintf(' external beta found ...'); end 
            IN.beta = InputParam.P{i}.BETAEXT; 
        else
            if isfield(InputParam.P{i},'SUBGROUP')
                IN.subgroup = InputParam.P{i}.SUBGROUP;
                if VERBOSE;fprintf(' subgroup beta computation ...'); end
            end
        end

    case 2
        % method 2: COMBAT HARMONIZATION.
        if VERBOSE
            fprintf('\tComBat harmonization ...');
            if IN.MCOVARUSE == 1; fprintf(' retaining variance from %d variables...', size(IN.TrMod, 2)); end
            if IN.MCOVARLABEL == 1; fprintf(' including label...'); end
        end
           
        if isfield(InputParam.P{i},'SUBGROUP')
            IN.subgroup = InputParam.P{i}.SUBGROUP;
            if VERBOSE;fprintf(' subgroup computation of metrics...'); end
        end

    case 3
        % METHOD 3: DISPARATE IMPACT REMOVER. 
        if VERBOSE
            fprintf('\tDisparate impact remover ...');
            if IN.DISTYPE == 1
                fprintf(' using median distribution ...');
            else
                fprintf(' using mean distribution ...');
            end
        end

        if isfield(InputParam.P{i},'SUBGROUP')
            IN.subgroup = InputParam.P{i}.SUBGROUP;
            if VERBOSE;fprintf('subgroup computation of distributions ...'); end
        end

        if isfield(actparam,'opt')
            InputParam.P{i}.LAMBDA = actparam.opt;
        end
end

if paramfl && tsfl && isfield(TrParami,'beta') 
    % Out-of-sample mode, used TsCovars and stored beta
    tsproc = true;
elseif trfl 
    if actparam.synthfl
        if iscell(SrcParam.covarsSyn)
            IN.TrCovars = [ InputParam.P{i}.TrCovars; SrcParam.covarsSyn{actparam.j}(:, IN.COVAR)]; 
        else
            IN.TrCovars = [ InputParam.P{i}.TrCovars; SrcParam.covarsSyn(:, IN.COVAR) ]; 
        end
    end
    [InputParam.Tr, TrParami] = nk_PartialCorrelationsObj(InputParam.Tr, IN);
    if tsfl, tsproc = true; end
else
    if VERBOSE;fprintf('not performed'); end
end

if tsproc
    if copy_ts
        InputParam.Ts{1} = InputParam.Tr;
        TrParami.copy_ts = true;
        Ts = nk_PartialCorrelationsObj(InputParam.Ts(2:end), TrParami); 
        InputParam.Ts(2:end) = Ts;
    else
        InputParam.Ts = nk_PartialCorrelationsObj(InputParam.Ts, TrParami); 
    end
end 

if isfield(CALIB,'calibflag') && CALIB.calibflag && isfield(InputParam, 'C')
    INcalib = IN;
    INcalib.TrCovars = CALIB.covars;
    INcalib = rmfield(INcalib, 'subgroup');
    %if isfield(InputParam, 'C') %CALIB.calibflag && CALIB.preprocstep > i
    InputParam.C = nk_PartialCorrelationsObj(InputParam.C, INcalib);
    %end
end

end

% =================================================================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_elimzero(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.PRUNE.perc = actparam.opt;
end

if paramfl && tsfl && isfield(TrParami,'NonPruneVec') 
    tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tAttribute pruning ...'); end
    [InputParam.Tr, TrParami] = nk_PerfElimZeroObj(InputParam.Tr, InputParam.P{i}.PRUNE);
    % All subsequent processing steps that use fixed column indices have to
    % be adjusted to the pruned matrix, except for the case of
    % dimensionality reduction preceeding these steps.
    InputParam = adjust_featuremasks(i, InputParam, TrParami.NonPruneVec);
    if tsfl, tsproc = true; end
end

if tsproc
    InputParam.Ts = nk_PerfElimZeroObj(InputParam.Ts, TrParami);
    InputParam.P{i,1}.INDPruneVec = find(~TrParami.NonPruneVec);
end
end

% =================================================================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_rankfeat(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

i  = actparam.i;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.RANK.opt = actparam.opt;
end
if isfield(InputParam,'Yw') && strcmp(InputParam.P{i}.RANK.algostr,'extern')
   InputParam.P{i}.RANK.EXTERN = InputParam.Yw;
end

RANK = InputParam.P{i}.RANK;

if ~isfield(TrParami,'W') || isempty(TrParami.W)
    % Modification allows to compute weights in a user-defined subgroup
    % (RANK.curglabel)
    if isfield(InputParam.P{i}.RANK,'curglabel')
        RANK.curlabel = InputParam.P{i}.RANK.curlabel(InputParam.P{i}.RANK.curglabel,actparam.curlabel);
        Y = InputParam.Tr(InputParam.P{i}.RANK.curglabel,:);
    else
        RANK.curlabel = InputParam.P{i}.RANK.curlabel(:,actparam.curlabel);
        Y = InputParam.Tr;
    end
    if actparam.synthfl
        % To do: synthetic data generation needs to create also
        % user-defined RANK.curlabel for artificial data. Currently only
        % synthetic target label can be created (TrainLabelSyn).
        if iscell(SrcParam.TrainLabelSyn)
            RANK.curlabel = [RANK.curlabel; SrcParam.TrainLabelSyn{actparam.j}(:,actparam.curlabel)]; 
        else
            RANK.curlabel = [RANK.curlabel; SrcParam.TrainLabelSyn(:,actparam.curlabel)]; 
        end
    end
    if VERBOSE;fprintf('\tCompute Feature Weighting ...'); end
    TrParami = nk_PerfFeatRankObj(Y, RANK);
    actparam.RANK.W{actparam.curlabel} = TrParami.W;
end

end

% =================================================================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_extfeat(SrcParam, InputParam, TrParam, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.W_ACT.opt = actparam.opt;
end

% Extract weight vector from last TrParam
if isfield(actparam,'j'); j = actparam.j; end
InputParam.P{i}.plsdev = false; 

if exist('j','var')
    InputParam.P{i}.W = TrParam{i-1}{j, actparam.curlabel}.W;
else
    if iscell(TrParam{i-1})
        InputParam.P{i}.W = TrParam{i-1}{actparam.curlabel}.W;
    else
        InputParam.P{i}.W = TrParam{i-1}.W;
    end
end

if VERBOSE,fprintf('\tRanking-based feature extraction ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    TrParami = InputParam.P{i};
    [ InputParam.Tr, TrParami ] = nk_PerfWActObj(InputParam.Tr, TrParami );
    if isfield(actparam,'PRUNE')
        ind = sum(TrParami.ind,2)>1;
        actparam.PRUNE.NonPruneVec = actparam.PRUNE.NonPruneVec(ind);
    end
    if tsfl, tsproc = true; end
    if isfield(InputParam,'Yw') && isnumeric(InputParam.Yw) && ~isempty(InputParam.Yw) && isfield(TrParami,'ind') 
        Yw = [];
        for i=1:size(TrParami.ind,2)
            Yw = [ Yw InputParam.Yw(TrParami.ind(:,i))];
        end
        InputParam.Yw = Yw;
    end
end

if tsproc, InputParam.Ts = nk_PerfWActObj(InputParam.Ts, TrParami ); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_remvarcomp(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.REMVARCOMP.corrthresh = nk_ReturnParam('CorrThresh', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 
    InputParam.P{i}.REMVARCOMP.DR.dims = nk_ReturnParam('Dims', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 
    %InputParam.P{i}.REMVARCOMP.DR.tolerance = nk_ReturnParam('fastICA_tol', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 

    %InputParam.P{i}.REMVARCOMP.DR.max_iters = nk_ReturnParam('fastICA_max_iter', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 
    
end

% We need the template parameters for the optional Procrustes rotation
if isfield(actparam,'Templ') && isfield(actparam.Templ,'mpp')
   InputParam.P{i}.REMVARCOMP.DR.mpp_template = actparam.Templ.mpp;
end

if VERBOSE,fprintf('\tExtract variance components ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    [ InputParam.Tr, TrParami ] = nk_PerfAdjForCovarsUsingPCAObj(InputParam.Tr, InputParam.P{i}.REMVARCOMP, InputParam.Tr);
    
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfAdjForCovarsUsingPCAObj(InputParam.Ts, TrParami ); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_devmap(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DEVMAP.cu = nk_ReturnParam('SPLS-cu', InputParam.P{i}.DEVMAP.Params_desc, actparam.opt); 
    InputParam.P{i}.DEVMAP.cv = nk_ReturnParam('SPLS-cv', InputParam.P{i}.DEVMAP.Params_desc, actparam.opt); 
end

%InputParam.P{i}.DEVMAP.loo = true;

if VERBOSE,fprintf('\tDeviation-based feature weighting ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl 
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    [ InputParam.Tr, TrParami ] = nk_PerfDevMapObj(InputParam.Tr, InputParam.P{i});
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfDevMapObj(InputParam.Ts, TrParami ); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_graphSparsity(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.GRAPHSPARSITY.p = actparam.opt;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tGraph sparsity thresholding ...'); end
    [InputParam.Tr, TrParami] = cv_graph_PerfSparsityThres(InputParam.Tr, InputParam.P{i}.GRAPHSPARSITY);

    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = cv_graph_PerfSparsityThres(InputParam.Ts, TrParami); end
end

function [SrcParam, InputParam, TrParami, actparam ] = act_graphMetrics(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.GRAPHMETRICS.p = actparam.opt;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tGraph metric computation ...'); end
    [InputParam.Tr, TrParami] = cv_graph_PerfGraphMetrics(InputParam.Tr, InputParam.P{i}.GRAPHMETRICS);
    % All 
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = cv_graph_PerfGraphMetrics(InputParam.Ts, TrParami); end
end

function [SrcParam, InputParam, TrParami, actparam ] = act_graphComputation(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE CALIB
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.GRAPHCONSTRUCTION.p = actparam.opt;
end

% if calibration data exists, replace refgroup with C
if CALIB.calibflag && isfield(InputParam, 'C') && strcmp(InputParam.P{i}.GRAPHCONSTRUCTION.method, "Normative network + 1") && InputParam.P{i}.GRAPHCONSTRUCTION.refGroup == -1
    InputParam.P{i}.GRAPHCONSTRUCTION.refGroup = InputParam.C;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tGraph computation ...'); end
    [InputParam.Tr, TrParami] = cv_graph_PerfGraphConstruction(InputParam.Tr, InputParam.P{i}.GRAPHCONSTRUCTION);
    % All 
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = cv_graph_PerfGraphConstruction(InputParam.Ts, TrParami); end
end
% =========================================================================

function [SrcParam, InputParam, TrParami, actparam ] = act_customPreproc(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.CUSTOMPREPROC.p = actparam.opt;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tCustom preprocessing function ...'); end
    [InputParam.Tr, TrParami] = cv_perfCustomPreproc(InputParam.Tr, InputParam.P{i}.CUSTOMPREPROC);

    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = cv_perfCustomPreproc(InputParam.Ts, TrParami); end
end
% =========================================================================

function [SrcParam, InputParam, TrParami, actparam ] = act_JuSpace(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.JUSPACE.p = actparam.opt;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tComputing neurotransmitter correlations ...'); end
    INDPruneVec = [];
    for prs = 1:size(InputParam.P,1)
        fields_curr = fieldnames(InputParam.P{prs,1});
        if any(contains(fields_curr,'INDPruneVec'))
            INDPruneVec = InputParam.P{prs,1}.INDPruneVec;
        end
    end
    [InputParam.Tr, TrParami] = cv_perfJuSpace(InputParam.Tr,InputParam.P{i}.JUSPACE,INDPruneVec);

    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = cv_perfJuSpace(InputParam.Ts, TrParami, INDPruneVec); end
end
% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_ROImeans(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE CALIB
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.ROIMEANS.p = actparam.opt;
end

if paramfl && tsfl 
     tsproc = true;
elseif trfl
    if VERBOSE;fprintf('\tROI means computation ...'); end
    [InputParam.Tr, TrParami] = cv_perfROImeans(InputParam.Tr, InputParam.P{i}.ROIMEANS);%,PREPROC,{InputParam.P{1:i-1}}); % prevPREPROC %prevP

    if tsfl, tsproc = true; end
end
% 
if tsproc, InputParam.Ts = cv_perfROImeans(InputParam.Ts, TrParami); end
if CALIB.calibflag && isfield(InputParam, 'C') %&& CALIB.preprocstep > i 
    InputParam.C = cv_perfROImeans(InputParam.C, InputParam.P{i}.ROIMEANS);
end
end
% =========================================================================
function [InputParam, SrcParam] = perform_adasyn(InputParam, SrcParam)
global SVM MODEFL

if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
    switch MODEFL
    case 'classification'
        TrL = SrcParam.oMultiTrainLabel;
    case 'regression'
        error('ADASYN works only for classification models');
    end
    if numel(unique(TrL))>2; error('ADASYN works currently only for binary classification'); end
    TrL(TrL>1)=0;
    try
        [InputParam.Tr, SrcParam.MultiTrainLabel] = nk_PerfADASYN( InputParam.Tr, TrL, SVM.ADASYN);  
    catch
        fprintf('error')
    end
else
    SrcParam.MultiTrainLabel = SrcParam.oMultiTrainLabel;
end

end

% =========================================================================
function InputParam = adjust_featuremasks(i, InputParam, NonPruneVec)

% Loop through the down-stream (i+z) preprocessing steps 
for z=i+1:numel(InputParam.P)
    switch InputParam.P{z}.cmd
        case 'reducedim'
            % If some sort of dimensionality reduction / factorization
            % method is performed we don't need to care anymore for
            % this problem
            break
        case 'impute' 
            if isfield(InputParam.P{z}.IMPUTE,'blockind') && ~isempty(InputParam.P{z}.IMPUTE.blockind)
                InputParam.P{z}.IMPUTE.blockind = InputParam.P{z}.IMPUTE.blockind(NonPruneVec);
            end
        case 'correctnuis'
            if isfield(InputParam.P{z},'BETAEXT') && ~isempty(InputParam.P{z}.BETAEXT)
                InputParam.P{z}.BETAEXT = InputParam.P{z}.BETAEXT(NonPruneVec);
            end
            % 24/01/2024: addition of featind correction
            if isfield(InputParam.P{z},'featind') && ~isempty(InputParam.P{z}.featinf)
                InputParam.P{z}.featind = InputParam.P{z}.featind(NonPruneVec);
            end
    end
end

end
