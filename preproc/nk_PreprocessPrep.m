function [act, analdim, p, GridAct, mapY, strout] = nk_PreprocessPrep( act, analdim, GridAct, p, parentstr)
% =========================================================================
% FORMAT function [act, analdim, p, GridAct, mapY, strout] = ...
%                   nk_PreprocessPrep( act, analdim, GridAct, p, parentstr)
% =========================================================================
% This function allows the interactive and batch use of the preprocessing
% module of NM.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 05/2024

global PREPROC MODEFL CV DR SAV RAND USEPARAMEXIST FUSION TEMPL CALIB MULTI STACKING NM OCTAVE CALIBUSE CVPOS
clc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('p','var') || isempty(p)
    p = struct('saveproc',      true, ...
                'ovrwrtfl',     2, ...
                'writefl',      2, ...
                'useparamfl',   2, ...
                'HideGridAct',  false, ...
                'batchflag',    false);
end

if p.batchflag == 2
    NM = p.NM;  
    CV = p.CV; 
end

yesnostr = {'yes','no'};

% Select analysis to work on

if ~isempty(analdim)
    if numel(analdim)<2
        AnalSelStr = sprintf('Analysis %g', analdim);
    else
        if ~p.HideGridAct, cvequalstr = 'same-size CV structures'; else, cvequalstr = 'different CV structures'; end
        AnalSelStr = sprintf('%g Analyses: %s [ %s ]',numel(analdim), strjoin(cellstr(num2str(analdim'))',', '), cvequalstr);
    end
else
    AnalSelStr = 'none selected';
end

if numel(NM.analysis) > 1 
    menustr = ['Select analysis to operate on [ ' AnalSelStr ' ]|']; menuact = 1;
else
    analdim = 1; menuact = []; menustr = [];
    analysis = NM.analysis{analdim}; nk_SetupGlobalVariables(analysis.params, 'setup_main', 0); 
    [operms,ofolds] = size(CV.TrainInd);
    if (~exist('GridAct','var') || isempty(GridAct)) || ~isequal(size(GridAct), [operms ofolds]) 
        GridAct = true(operms,ofolds);
    end
end

% Overwriting PreprocData mats?
menustr = [menustr ...
            'Overwrite existing preprocessed data files [ ' yesnostr{p.ovrwrtfl} ' ]|' ...
            'Use existing preprocessing parameter files, if available [ ' yesnostr{p.useparamfl} ' ]'];
menuact = [menuact 2 3];

% Writing parameters to disk?
if p.useparamfl == 2
    menustr = [menustr '|Write out computed parameters to disk (may require A LOT of disk space) [ ' yesnostr{p.writefl} ' ]'];
    menuact = [menuact 4];
end

% Select CV2 partitions
if ~isempty(analdim) && ~p.HideGridAct
    menustr = [menustr '|Select CV2 partitions to operate on [ ' num2str(sum(GridAct(:))) ' selected ] ']; menuact = [menuact 5];
end

% Allow user to proceed with the analysis
if (sum(GridAct(:)) || p.HideGridAct) && ~isempty(analdim)
    menustr = [menustr '|PROCEED >>>']; menuact = [menuact 6];
end

nk_PrintLogo
mestr = 'Preprocessing module run-time configuration'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
if ~p.batchflag, act = nk_input(mestr, 0, 'mq', menustr, menuact); end

switch act
    case 1
        showmodalvec=[]; tact = 1; tanaldim = analdim; brief=1;
        while tact>0, [ tact, tanaldim,~, showmodalvec, brief ] = nk_SelectAnalysis(NM, 0, 'MAIN INTERFACE >> PREPROCESSING MODULE', tanaldim, [], 0, showmodalvec, brief); end;
        if ~isempty(tanaldim) 
            analdim = tanaldim;
            analysis = NM.analysis{analdim}; nk_SetupGlobalVariables(analysis.params, 'setup_main', 0); 

            nA = numel(analdim);
            if nA>1
                AS = nk_GetAnalysisStatus(NM, analdim);
                if ~AS.betweenfoldpermequal_cv
                    p.HideGridAct = true; 
                else
                    [ operms, ofolds ] = size(NM.analysis{analdim(1)}.params.cv.TrainInd);
                end
            else
                p.HideGridAct = false;
                [ operms, ofolds ] = size(NM.analysis{analdim}.params.cv.TrainInd);
            end
            
            if ~p.HideGridAct && (~exist('GridAct','var') || isempty(GridAct) || ~isequal(size(GridAct), [operms ofolds])) 
                GridAct = false(operms,ofolds);
            end
        end
    case 2
        if p.ovrwrtfl == 1, p.ovrwrtfl = 2; elseif p.ovrwrtfl == 2, p.ovrwrtfl = 1; end
    case 3
        if p.useparamfl == 1, p.useparamfl = 2; elseif p.useparamfl == 2, p.useparamfl = 1; end
    case 4
        if p.writefl == 1, p.writefl=2; elseif p.writefl == 2, p.writefl = 1; end
    case 5
        %% Get info from user which CV2 partitions to operate upon?   
        [operms,ofolds] = size(NM.analysis{analdim(1)}.params.cv.TrainInd);
        t_act = 1; while t_act > 0 && t_act < 10, [ t_act, GridAct ] = nk_CVGridSelector(operms, ofolds, GridAct, 0); end

    case {6,7}
        nA=numel(analdim);

        for i=1:nA

            analysis = NM.analysis{analdim(i)};
            clc
            fprintf('**********************************\n')
            fprintf('**  PREPROCESSING OPTIMIZATION  **\n')
            fprintf('**********************************\n')
            fprintf('\nWorking on analysis #%g: %s', analdim(i), analysis.id);

            nk_SetupGlobalVariables(analysis.params, 'setup_main', 0); 
            if p.HideGridAct
                [ operms, ofolds] = size(NM.analysis{analdim(i)}.params.cv.TrainInd);
                GridAct = true(operms,ofolds);
            else
                [operms,ofolds] = size(CV.TrainInd);
            end

            NM.runtime.curanal = analdim(i);
            tdir = pwd; if isfield(analysis,'rootdir') && exist(analysis.rootdir,'dir'), tdir = analysis.rootdir; end
            if ~exist(analysis.rootdir,'dir'), mkdir(analysis.rootdir);end
            % Define modality-independent parameters of current analysis as global variables
            
            if ~isempty(STACKING) && STACKING.flag ==1
                act = 7; % Generate features for stacker
                F = 1; nM = 1;
            else
                if ~isempty(FUSION)        
                    F = analysis.params.TrainParam.FUSION.M;
                    nM = numel(F);
                    if FUSION.flag == 1, nM = 1; end
                else
                    F = 1; nM = 1;
                end
            end     
            
            matname = SAV.matname;
            if p.ovrwrtfl == 1
                OVRWRT = true;
            else
                OVRWRT = false;
            end
            
            %% Set further important program flags
            if  p.writefl == 1
                paramfl.write = 1; WRITEPARAM = true; 
            else
                paramfl.write = false; WRITEPARAM = false;
            end
    
            %% Set further important program flags
            if  p.useparamfl == 1
                paramfl.use_exist = true; USEPARAMEXIST = true;
            else
                paramfl.use_exist = false; USEPARAMEXIST = false;
            end
            
            CALIBUSE = false;
            kbin = 1;
            %% Define program parameters
            if strcmp(MODEFL,'classification') && RAND.Decompose ~= 9, kbin = length(CV.class{1,1}); end
            datid       = NM.id;
        
            switch act
                case 6

                 for j = 1:nM

                    paramfl.currmodal = j;
                    
                    %% Get Training / CV data (Y) & Build modality suffix
                    inp = nk_DefineFusionModeParams(NM, analysis, F, nM, j);
    
                    %fprintf('\nWORKING ON MODALITY #%g', j);
                    if numel(inp.PREPROC)>1
                        Y = inp.X(j).Y; 
                        if isfield(inp.X(j),'altY')
                            inp.AltY = inp.X(j).altY;
                        end
                        PREPROC = inp.PREPROC{j};
                        if isfield(inp.X(j),'Yw'); inp.iYw = inp.X(j).Yw; end
                    else
                        Y = inp.X.Y; 
                        if isfield(inp.X,'altY')
                            inp.AltY = inp.X(j).altY;
                        end
                        if isfield(inp.X,'Yw'); inp.iYw = inp.X.Yw; end
                        PREPROC = inp.PREPROC; 
                    end
    
                    DR = get_preproc_params(PREPROC);
    
                    % ======================== PRE-CV PROCESSING STEPS =======================
                    % No cross-validation is needed for these processing steps as these steps 
                    % are carried out separately for each subject, WITHOUT taking into account 
                    % information from other subjects
    
                    % LABEL SCALING (regression only)
                    labels = nk_LabelTransform(PREPROC, MODEFL, inp.labels);
    
                    % For imaging data: DATA FILTERING / SMOOTHING / RESLICING (in the future)
                    Y = nk_PerfSpatFilt( Y, PREPROC, inp.X );
                    inp = nk_SmoothMask( PREPROC, inp );

                    CALIB.calibflag = false;
                    % Check whether calibration data is available 
                    if isfield(inp, 'C') && inp.C{1,1}.calibflag
                        CALIB.calibflag = inp.C{1,1}.calibflag;
                        P = inp.X;
                        CYfile = inp.C{1,1}.Y; 
                        load(CYfile, 'CY');
                        inp.C{1,1}.Y = CY;
                        inp.C{1,1}.Y = nk_PerfSpatFilt(inp.C{1,1}.Y, PREPROC, P);
                    end

                    % ======================== PREPROCESSING PIPELINE ========================
                    % These stepps require cross-validation as they require group-level
                    % information flows
                    for ix=1:operms % Loop through CV2 perms
    
                        for jx=1:ofolds % Loop through CV2 folds
    
                            if ~GridAct(ix,jx), continue; end
    
                            inp.f = ix; inp.d = jx; inp.nclass = kbin;
                            CVPOS.CV2p = ix; CVPOS.CV2f = jx;
							
                            %% Construct output filename(s) and check for their existence
                            strout = nk_Preprocess_StrCfg(ix, jx);
                            % Run preprocessing on the data matrix
                           
                            savnamY = [matname strout inp.varstr '_PreprocData_ID' datid];
                            savmatY  = fullfile(tdir,[savnamY '.mat']);
                            savnamP = [matname strout inp.varstr '_PreprocDataParam_ID' datid];
                            savmatP = fullfile(tdir,[savnamP '.mat']);
                            paramfl.pth = savmatP; 

                            if ~OVRWRT
                                flg=0;
                                if exist(savmatY,'file'), fprintf('\n Training / CV file detected:\n%s\nDo not overwrite.', savnamY); flg=1; end
                                if WRITEPARAM, if exist(savmatP,'file'), fprintf('\n Parameter file detected:\n%s\n%s\nDo not overwrite.', savnamP); end; end
                                if flg, continue; end
                            end

                            inp.maindir = analysis.rootdir;
                            inp.id = datid;

                            fprintf('\n\n'); fprintf('********************** CV2 partition [%g, %g] ********************** ',ix,jx)
    
                            %% Prepare parameter container, if paramfl = true
                            paramfl.found = 0;
                            %paramfl.writeCV1 = true;
                            if paramfl.use_exist
                                try 
                                    load(paramfl.pth);
                                    paramfl.Param = Param; clear Param;
                                    paramfl.found = 1;
                                catch
                                    paramfl.found = 0;
                                end
                            end
                            
                            %% For factorization methods: TEMPLATE MAPPING   
                            if isfield(PREPROC,'TEMPLPROC') && ~isempty(PREPROC.TEMPLPROC) && PREPROC.TEMPLPROC
                                paramfl.templateflag = true;
                            end
                            
                            
                            %% Perform preprocessing
                            [mapY, Param] = nk_PerfPreprocess(Y, inp, labels, paramfl);
    
                            outfold = jx; outperm = ix;
                            
                            if p.saveproc
                                fprintf('\n* Saving TRAINING / CV data:\n%s\n', savnamY)
                                if OCTAVE
                                    save(savmatY,'mapY','datid', 'outperm', 'outfold')
                                else
                                    save(savmatY,'mapY','datid', 'outperm', 'outfold', '-v7.3')
                                end
                                if WRITEPARAM
                                    fprintf('\n* Saving Parameter data:\n%s\n', savnamP)
                                    if OCTAVE
                                        save(savmatP, 'Param','datid', 'outperm', 'outfold')
                                    else
                                        save(savmatP, 'Param','datid', 'outperm', 'outfold', '-v7.3')
                                    end
                                end
                                clear FEAT Param mapY id dimension savmat savnam strout volnum sigflag clustflag signum clustnum
                            end
                            if exist('TEMPL','var') && ~isempty(TEMPL), clear TEMPL; end
                        end
                    end
                end
    
            case 7

                nk_SetupGlobalVariables(analysis.params, 'setup_strat', 1);                 

                %% Get Training / CV data (Y) & Build modality suffix
                inp = nk_DefineFusionModeParams(NM, analysis, F, nM, 1);
                inp.nM = nM;
                inp.multiflag = 0; if ~isempty(MULTI) && MULTI.flag, inp.multiflag = nk_input('Select models at multi-group optimum', 0,'yes|no',[1,0],1); end
                if isfield(NM,'covars'); inp.covars = NM.covars; else, inp.covars = []; end
                inp.covars_oocv = []; inp.ll=1;
                
                % LABEL SCALING (regression only)
                labels = nk_LabelTransform(PREPROC, MODEFL, inp.labels);
                
                for ix=1:operms % Loop through CV2 perms
    
                        for jx=1:ofolds % Loop through CV2 folds
    
                            if ~GridAct(ix,jx), continue; end
                            outfold = jx; outperm = ix;
                            inp.f = ix; inp.d = jx; inp.nclass = kbin;
                            %% Construct output filename(s) and check for their existence
                            strout = nk_Preprocess_StrCfg(ix, jx);
                            savnamY = [matname strout '_PreprocDataMeta_ID' NM.id];
                            savmatY  = fullfile(tdir,[savnamY '.mat']);
                            savnamP = [matname strout inp.varstr '_PreprocDataParamMeta_ID' datid];
                            savmatP = fullfile(tdir,[savnamP '.mat']);
                          
                            paramfl.found = 0;
                            if paramfl.use_exist
                                try 
                                    load(paramfl.pth);
                                    paramfl.Param = Param; clear Param;
                                    paramfl.found = 1;
                                catch
                                    paramfl.found = 0;
                                end
                            end
                            for curlabel = 1:size(labels,2)
                                inp.curlabel=curlabel;
                                mapY(curlabel) = nk_PerfPreprocessMeta(inp, labels(:,curlabel), paramfl);
                            end                                 
                            outfold = jx; outperm = ix;
                            
                            if p.saveproc
                                fprintf('\n* Saving META-LEVEL TRAINING / CV data:\n%s\n', savnamY)
                                if OCTAVE
                                    save(savmatY,'mapY','datid', 'outperm', 'outfold')
                                else
                                    save(savmatY,'mapY','datid', 'outperm', 'outfold', '-v7.3')
                                end
                                if WRITEPARAM
                                    fprintf('\n* Saving Parameter data:\n%s\n', savnamP)
                                    if OCTAVE
                                        save(savmatP, 'Param','datid', 'outperm', 'outfold')
                                    else
                                        save(savmatP, 'Param','datid', 'outperm', 'outfold', '-v7.3')
                                    end
                                end
                                clear FEAT Param mapY id dimension savmat savnam strout volnum sigflag clustflag signum clustnum
                            end
                            inp.ll = inp.ll + 1;
                        end
                end
            end
            nk_SetupGlobalVariables(NM.analysis{analdim(i)}.params, 'clear', 0);
        end
end

%fprintf('\n* Preprocessing finished.')

function [DR, BINMOD, FEATSEL, CLUST, COVAR] = get_preproc_params(PREPROC)

DR=[]; FEATSEL=[]; CLUST=[]; COVAR=[]; BINMOD=[]; 
if isfield(PREPROC,'FEATSEL'), FEATSEL = PREPROC.FEATSEL; end
if isfield(PREPROC,'BINMOD'), BINMOD = PREPROC.BINMOD; end
if isfield(PREPROC,'ACTPARAM')
    for i=1:numel(PREPROC.ACTPARAM)
       if isfield(PREPROC.ACTPARAM{i},'DR')
           DR = PREPROC.ACTPARAM{i}.DR;
       elseif isfield(PREPROC.ACTPARAM{i},'CLUST')
           CLUST = PREPROC.ACTPARAM{i}.CLUST;
       elseif isfield(PREPROC.ACTPARAM{i},'COVAR')
           COVAR = PREPROC.ACTPARAM{i}.COVAR;  
       end
    end
end
