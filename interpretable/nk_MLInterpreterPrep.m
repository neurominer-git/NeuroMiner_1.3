function [ act, NM, inp ] = nk_MLInterpreterPrep(NM, act, inp, parentstr)
% =========================================================================
% FORMAT [act, inp] = nk_MLInterpreterPrep(act, inp, parentstr)
% =========================================================================
% Runtime module of interpretable ML module
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, last modified 09/2022

global CV 
 
OverWriteStr = []; GridSelectStr = []; LoadModelsStr = []; LoadParamsStr = []; LoadStr = []; SaveStr = []; SaveCV1Str = [];
OverWriteAct = []; GridSelectAct = []; LoadModelsAct = []; LoadParamsAct = []; LoadAct = []; SaveAct = []; SaveCV1Act = [];
DATASCRAM = false; if isfield(NM.defs,'data_scrambled') && ~isempty(NM.defs.data_scrambled), DATASCRAM = NM.defs.data_scrambled; end

as = nk_GetAnalysisStatus(NM); complvec = find(as.completed_analyses);

% Initialize runtime parameters 
if ~exist('inp','var') || isempty(inp)
    if isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked && isfield(NM,'OOCV')
        oocvflag = true;
        oocvind = 1;
        OO = NM.OOCV{1};
    else
        oocvflag = false;
        oocvind = [];
        OO = [];
    end

    inp = struct(  'analind', complvec(1), ...  % Index to analysis
                    'oocvflag', oocvflag, ...   % Are we in OOCV mode?
                    'oocvind', oocvind, ...     % Index to OOCV data container
                    'OO', OO, ...               % User-defined OOCV data container (start with oocvind = 1)
                    'lfl', 1, ...               % 1 = compute from scratch | 
                    ...                         % 2 = use existing (allowing the user to specify OOCVdatamats)
                    'ovrwrt', 2, ...            % if lfl == 1 ==> 1 = overwrite existing OOCVdatamats 
                    ...                         % 2 = do not overwrite (use existing OOCVdatamats automatically)
                    'saveparam', 2, ...         % if loadparam == 2=> 1 = save OOCV processing parameters (preprocessing / models)
                    'ovrwrtperm', 2, ...        % Overwrite permutation file.
                    'recompute_estimates', 2,...% Recompute prediction change estimates
                    ...                         % 2 = do not save parameters to disk
                    'saveCV1', 2, ...           % if loadparam == 2 && saveparam ==1 => 1 = save large OOCV processing containers at the CV1 level 
                    ...                         % 2 = operate at CV2 level
                    'loadparam', 2, ...         % 1 = load existing optpreproc and/or optmodel parameters from disk
                    ...                         % 2 = recompute parameters
                    'selectsamples', 1, ...     % 1 = All
                    ...                         % 2 = User-defined
                    'sampleID', {{}}, ...         % cell array of strings defining IDs for which to run MLI if user-defined 
                    'HideGridAct', false, ...
                    'batchflag', 0);            % 1 = Run in batchmode (without graphics outputs)
                                                % 0 = run in interactive mode
                    % 'sampleID', {}, ...         % cell array of strings defining IDs for which to run MLI if user-defined 


end
% Prepare variables for menu creation
na_str = '?'; inp.datatype = 'MLIdatamat';

%% Configure menu
% Select analysis
if numel(NM.analysis)>1
    if numel(inp.analind)<2
        AnalSelStr = sprintf('Analysis %g', inp.analind);
    else
        if ~inp.HideGridAct, cvequalstr = 'same-size CV structures'; else, cvequalstr = 'different CV structures'; end
        AnalSelStr = sprintf('%g Analyses: %s [ %s ]',numel(inp.analind), strjoin(cellstr(num2str(inp.analind'))',', '), cvequalstr);
    end 
    AnalSelectStr = sprintf('Choose analysis to work on [ %s ]|', AnalSelStr);    AnalSelectAct=1;    
else
    AnalSelectStr = ''; AnalSelectAct = [];
end
analysis      = NM.analysis{inp.analind(1)}; 

if isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked
    % Select independent test data container
    inp.oocvflag = true;
    if isfield(inp,'oocvind') 
        OOCVSelStr = sprintf('New data #%g: %s', inp.oocvind, inp.OO.desc); 
    else 
        OOCVSelStr = na_str; 
    end
    OOCVSelectStr = sprintf('Choose independent data to work on [ %s ]|', OOCVSelStr);                                    OOCVSelectAct = 2;   
    if DATASCRAM, inp.loadparam = 1; inp.saveparam = 2; end
else
    % Otherwise interpret CV2 validation data
    inp.oocvflag = false;
    OOCVSelectStr = '';
                                                                                                                          OOCVSelectAct = [];
end

if ~isempty(analysis)
    
    % Initialize global parameters for the selected analysis
    nk_SetupGlobalVariables(analysis.params, 'setup_main', 0); 
    
    % Compute from scratch or use pre-computed datamats ?
    LFL_opts        = {'Compute from scratch',sprintf('Use precomputed %s',inp.datatype)};                                      
    ModeStr         = sprintf('Operation mode of ML interpreter module [ %s ]|',LFL_opts{inp.lfl});                       ModeAct = 3;
    OVRWRT_opts     = {'Overwrite existing','Do not overwrite'}; 
    if inp.lfl == 1
        % from scratch      
        OverWriteStr = sprintf('Overwrite existing %s files [ %s ]|', inp.datatype, OVRWRT_opts{inp.ovrwrt}) ;            OverWriteAct = 4; 
    else
        % precomputed
        nMatFiles = na_str; 
        if isfield(inp,'matfiles') && ~isempty(inp.matfiles) 
            selGrid = ~cellfun(@isempty,inp.matfiles); inp.GridAct = selGrid;
            nMatFiles = sprintf('%g selected', sum(selGrid(:))); 
        end     
        OverWriteStr = sprintf('Specify %s files [ %s ]|', inp.datatype, nMatFiles);                                      OverWriteAct = 4; 
    end
    
    % Overwrite permutation file
    OverWritePermStr = sprintf('Overwrite existing permutations file [ %s ]|', OVRWRT_opts{inp.ovrwrtperm}) ;             OverWritePermAct = 5; 
    
    % Recompute prediction change estimate
    LOAD_opts        = {'yes', 'no'}; 
    RecomputeChangeStr = sprintf('Recompute prediction change estimates [ %s ]|', LOAD_opts{inp.recompute_estimates}) ;   RecomputeChangeAct = 6; 

    % Retrieve CV2 partitions to operate on
    if ~isfield(inp,'GridAct'), inp.GridAct = analysis.GDdims{1}.GridAct; end
    if ~inp.HideGridAct
        GridSelectStr = sprintf('Select CV2 partitions to operate on [ %g selected ]|',  sum(inp.GridAct(:)));                  GridSelectAct = 7;
    else
        GridSelectStr =''; GridSelectAct=[];
    end
    
    % Configure loading of pre-existing parameters and models
    if inp.saveparam == 2 && inp.lfl == 1
        
        if ~DATASCRAM
            LoadStr = sprintf('Use saved pre-processing params and models [ %s ]|', LOAD_opts{inp.loadparam});            LoadAct = 8;
        end
        if inp.loadparam == 1
            if isfield(inp,'optpreprocmat') 
                selGridPreproc = ~cellfun(@isempty,inp.optpreprocmat);
                nParamFiles = sprintf('%g files selected', sum(selGridPreproc(:))); 
            else
                nParamFiles = na_str; 
            end
            LoadParamsStr = sprintf('Select preprocessing parameter files [ %s ]|' ,nParamFiles);                         LoadParamsAct = 10;
            if isfield(inp,'optmodelmat') 
                selGridModel = ~cellfun(@isempty,inp.optmodelmat);
                nModelFiles = sprintf('%g files selected', sum(selGridModel(:))); 
            else 
                nModelFiles = na_str; 
            end
            LoadModelsStr = sprintf('Select model files [ %s ]|',nModelFiles);                                            LoadModelsAct = 11;
        end
    end
    
    % If loading of pre-existing models and params is not chosen, allow to
    % save the computed params and models to disk
    if inp.loadparam == 2 && inp.lfl == 1
        SAVE_opts       = {'yes', 'no'};   
        SaveStr = sprintf('Save pre-processing params and models to disk [ %s ]|', SAVE_opts{inp.saveparam});             SaveAct = 9;
    end

    % Select sample to run MLI for
    SAMPLE_opts        = {'All', 'User-defined'}; 
    SelectSamplesStr   = sprintf('Define for which sample(s) to interpret the model [ %s ]|', SAMPLE_opts{inp.selectsamples}) ;   
    SelectSampleAct    = 14; 

    % Select sample(s) to run MLI for
    if isequal(SAMPLE_opts{inp.selectsamples}, 'User-defined')
        if isempty(inp.sampleID)
            sampleIDstr = '';
        elseif numel(inp.sampleID) > 4
            sampleIDstr = [num2str(numel(inp.sampleID)) ' IDs'];
        else
            sampleIDstr = strjoin(inp.sampleID);
        end
        SampleIDStr         = sprintf('Specify sample ID(s) for which to interpret the model [ %s ]|', sampleIDstr) ;   
    else
        SampleIDStr         = sprintf('');
    end
    SelectedSampleAct   = 15; 
end
 
 %% Build interactive menu
menustr = [ AnalSelectStr ...
            OOCVSelectStr ...
            ModeStr ...
            OverWriteStr ...
            OverWritePermStr ...
            RecomputeChangeStr ...
            GridSelectStr ...
            SaveStr ...
            SaveCV1Str ...
            LoadStr ...
            LoadParamsStr ... 
            LoadModelsStr ...
            SelectSamplesStr ...
            SampleIDStr];

menuact = [ AnalSelectAct ...
            OOCVSelectAct ...
            ModeAct ...
            OverWriteAct ...
            OverWritePermAct ...
            RecomputeChangeAct ...
            GridSelectAct ...
            SaveAct ...
            SaveCV1Act ...
            LoadAct ...
            LoadParamsAct ...
            LoadModelsAct ...
            SelectSampleAct];  

if ~isempty(SampleIDStr)
    menuact = [menuact SelectedSampleAct];
end

disallow = false;

%% Check whether all parameters are available
if (~sum(inp.GridAct(:)) && ~inp.HideGridAct) || isempty(inp.analind), disallow = true; end

if inp.loadparam == 1
    if ~isfield(inp,'optpreprocmat') || isempty(inp.optpreprocmat), disallow = true; end
    if ~isfield(inp,'optmodelmat') || isempty(inp.optmodelmat), disallow = true; end
end

if ~disallow, menustr = [menustr '|PROCEED >>>']; menuact = [menuact 12]; end 

%% Display menu and act on user selections
nk_PrintLogo
mestr = 'Model interpreter module run-time configuration'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
if ~inp.batchflag && act<16, act = nk_input(mestr, 0, 'mq', menustr, menuact); end

switch act
    
    case 0
        return
    % Select analysis
    case 1 
        showmodalvec = []; analind = inp.analind; 
        if length(NM.analysis)>1, t_act = 1; brief = 1;
            while t_act>0
                [t_act, analind, ~, showmodalvec , brief] = nk_SelectAnalysis(NM, 0, navistr, analind, [], 1, showmodalvec, brief, 7); 
            end
            if ~isempty(analind)
                %inp.analind = indanal(analind) ; 
                inp.analind = analind; 
            end
            nA = numel(inp.analind);
            if nA>1
                AS = nk_GetAnalysisStatus(NM, inp.analind);
                if ~AS.betweenfoldpermequal_cv
                    inp.HideGridAct = true; 
                else
                    inp.GridAct = NM.analysis{inp.analind(1)}.GDdims{1}.GridAct;
                    inp.HideGridAct = false;
                end
            else
                inp.HideGridAct = false;
                inp.GridAct = NM.analysis{inp.analind}.GDdims{1}.GridAct;
            end
        end
        
    % Select OOCV data if in OOCV mode
    case 2    
        [ NM, OO, oocvind ] = nk_SelectOOCVdata(NM, 1, 0);  
        if ~isempty(oocvind), inp.OO = OO; inp.oocvind = oocvind; end 
    case 3
         lfl = nk_input('Define run-time mode of ML interpreter module',0,'mq',strjoin(LFL_opts, '|'),[1,2],inp.lfl);
         if lfl, inp.lfl = lfl; end 
    % Overwrite?
    case 4
        switch inp.lfl
            case 1
                if inp.ovrwrt == 1, inp.ovrwrt = 2; elseif inp.ovrwrt  == 2, inp.ovrwrt = 1; end
            case 2
                if inp.oocvflag
                    tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
                else
                    tdir = create_defpath(NM.analysis{inp.analind});
                end
                inp.matfiles = nk_GenDataMaster(NM.id, inp.datatype, CV, [], tdir);
        end
    case 5
        if inp.ovrwrtperm == 1, inp.ovrwrtperm = 2; elseif inp.ovrwrtperm == 2, inp.ovrwrtperm = 1; end
    case 6
        if inp.recompute_estimates == 1, inp.recompute_estimates = 2; elseif inp.recompute_estimates == 2, inp.recompute_estimates = 1; end
    case 7
        [operms,ofolds] = size(CV.TrainInd);
        tact = 1; while tact > 0 && tact < 10, [ tact, inp.GridAct ] = nk_CVGridSelector(operms, ofolds, inp.GridAct, 0); end
    case 8
        if inp.saveparam == 1, inp.saveparam = 2; elseif inp.saveparam == 2,  inp.saveparam = 1; end
    case 9
        if inp.loadparam == 1, inp.loadparam = 2; elseif inp.loadparam == 2,  inp.loadparam = 1; end
    case 10
        tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
        optpreprocmat = nk_GenDataMaster(NM.id, 'OptPreprocParam', CV, [], tdir);
        if ~isempty(optpreprocmat), inp.optpreprocmat = optpreprocmat; end
    case 11
        tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
        optmodelmat = nk_GenDataMaster(NM.id, 'OptModel', CV, [], tdir);
        if ~isempty(optmodelmat), inp.optmodelmat = optmodelmat; end
    case {12,13}
        if inp.oocvflag
            inp.oocvname = sprintf('OOCV_%g',inp.oocvind);
        end
        nA = 1; if numel(inp.analind)>1, nA = numel(inp.analind); end
        for i=1:nA
            NM.runtime.curanal = inp.analind(i);
            inp = nk_GetAnalModalInfo_config(NM, inp);
            if inp.HideGridAct, [ ix, jx ] = size(NM.analysis{inp.analind(i)}.params.cv.TrainInd); inp.GridAct = true(ix,jx); end
            inp.analysis_id = NM.analysis{inp.analind(i)}.id;
            inp.saveoptdir = [ NM.analysis{inp.analind(i)}.rootdir filesep 'opt' ];
            if inp.oocvflag
                NM.analysis{inp.analind(i)}.OOCV{inp.oocvind}.MLI = MLInterpreterPrep(NM, inp, NM.analysis{inp.analind(i)});
            else
                NM.analysis{inp.analind(i)}.MLI = MLInterpreterPrep(NM, inp, NM.analysis{inp.analind(i)});
            end
            nk_SetupGlobalVariables(NM.analysis{inp.analind(i)}.params, 'clear', 0); 
        end
        NM = rmfield(NM,'runtime'); 
    case 14
        if inp.selectsamples == 1, inp.selectsamples = 2; elseif inp.selectsamples == 2,  inp.selectsamples = 1; end
    case 15
        try
            inp.sampleID = strjoin(inp.sampleID);
        catch
            inp.sampleID = inp.sampleID;
        end
        while ~iscellstr(inp.sampleID)
        inp.sampleID = nk_input('Enter sample IDs as a cell array of strings (e.g. {''ID1'',''ID2''})', 0, 'e', inp.sampleID);
        end
        %make sure the cell array of strings is in one column
        inp.sampleID = reshape(inp.sampleID,size(inp.sampleID,1)*size(inp.sampleID,2),1);
end

function tdir = create_defpath(analysis, oocvind)
 
rootdir = analysis.GDdims{1}.RootPath;
if exist("oocvind","var") 
    if isfield(analysis,'OOCV') && ...
        numel(analysis.OOCV) >= oocvind && ...
        isfield(analysis.OOCV{oocvind},'RootPath')
        if iscell(analysis.OOCV{oocvind}.RootPath)
            tdir = analysis.OOCV{oocvind}.RootPath{1};
        else
            tdir = analysis.OOCV{oocvind}.RootPath;
        end
    else
        oocvdir = sprintf('OOCV_MLI_%g', oocvind);
        tdir = fullfile(rootdir, oocvdir);
    end
else
     tdir = fullfile(rootdir, 'MLI');
end

%
% =========================================================================
function MLIres = MLInterpreterPrep(dat, inp1, analysis)
global SAV MODEFL CV FUSION MULTILABEL MLI
if inp1.saveparam   == 2, inp1.saveparam    = 0; end
if inp1.loadparam   == 2, inp1.loadparam    = 0; end
if inp1.ovrwrt      == 2, inp1.ovrwrt       = 0; end
if inp1.lfl         == 1, inp1.analmode     = 0; else, inp1.analmode = 1; end
inp1.multiflag = 0;
F = 1; nF = 1;
if ~isempty(FUSION)        
    F = analysis.params.TrainParam.FUSION.M;
    nF = numel(F); if FUSION.flag < 3, nF = 1; end
    inp1.nF = nF;
end

if strcmp(MODEFL,'classification')
    inp1.nclass = length(CV.class{1,1});
else
    inp1.nclass = 1;
end

if inp1.oocvflag
    if isfield(inp1.OO,'label') && ~isempty(inp1.OO.label)
        inp1.LabelCV     = dat.label; 
        inp1.labelOOCV   = inp1.OO.label; 
    end
    inp1.cases_oocv      = inp1.OO.cases;
    inp1.nOOCVsubj       = numel(inp1.OO.cases);
end
inp1.cases           = dat.cases;
inp1.id              = dat.id;
stranalysis          = SAV.matname;
inp1.ngroups         = numel(unique(dat.label(~isnan(dat.label))));

if isfield(inp1.OO,'groups') && numel(inp1.OO.groups)==numel(inp1.OO.cases)
    inp1.groupind = inp1.OO.groups;
    if isfield(inp1.OO,'groupnames')
        inp1.groupnames = inp1.OO.groupnames;
    end
end

if inp1.oocvflag
    if isfield(inp1,'targdir') 
        inp1.rootdir = fullfile(inp1.targdir, [ inp1.oocvname '_MLI'] );
    elseif isfield(analysis,'rootdir') && exist(analysis.rootdir,'dir')
        inp1.rootdir = fullfile(analysis.rootdir,analysis.params.TrainParam.SVM.prog, [ inp1.oocvname '_MLI']);
    else
        inp1.rootdir = fullfile(pwd,analysis.params.TrainParam.SVM.prog, [ inp1.oocvname '_MLI']);
    end
else
    inp1.rootdir = fullfile(analysis.rootdir,analysis.params.TrainParam.SVM.prog,'MLI');
end
inp1.maindir = analysis.rootdir;
if ~exist(inp1.rootdir,'dir'), mkdir(inp1.rootdir); end
nl = nk_GetLabelDim(MULTILABEL);
inp1.MLI = MLI;

% Loop through modalities
for i = 1:inp1.nF
    
    % **************************** ANALYSIS SETUP *****************************
    inp2 = nk_DefineFusionModeParams(dat, analysis, F, nF, i, inp1.oocvind);
    inp2.labels = analysis.params.label.label;
    inp = catstruct(inp1,inp2);
    inp.MLI.Modality = MLI.Modality(inp.tF);
    inp.loadGD = true;

    for j = 1:nl
	    
        if inp.oocvflag
            Resultsfile = fullfile(inp.rootdir,[stranalysis inp.varstr '_t' num2str(j) '_OOCV-' num2str(inp.oocvind) '-MLIresults_ID' dat.id '.mat']);
        else
            Resultsfile = fullfile(inp.rootdir,[stranalysis inp.varstr '_t' num2str(j) '_MLIresults_ID' dat.id '.mat']);
        end
        inp.multlabelstr = '';  if MULTILABEL.flag, inp.multlabelstr = sprintf('_t%g',j); end
        if exist(Resultsfile,'file') && inp.ovrwrt==2
	        fprintf('\nLoading:\n%s',Resultsfile);
	        load(Resultsfile)
        else
	        if MULTILABEL.flag && MULTILABEL.dim>1
		        fprintf('\n\n');fprintf('====== Working on label #%g ====== ',j);
		        inp.curlabel = j;
	        else
		        inp.curlabel = 1;
	        end
	        if strcmp(MODEFL,'classification')
		        [ijMLI.BinResults, ijMLI.FileNames, ijMLI.RootPath] = nk_MLInterpreter(inp);
			else
                [ijMLI.RegrResults, ijMLI.FileNames, ijMLI.RootPath] = nk_MLInterpreter(inp);
	        end
	        fprintf('\nSaving:\n%s',Resultsfile);
	        save(Resultsfile,'ijMLI', 'MLI');	
        end
        
        switch MODEFL
            case 'classification'
                MLIres.Label(i,j) = ijMLI.BinResults;  
            case 'regression'
                MLIres.Label(i,j) = ijMLI.RegrResults;
        end
        MLIres.FileNames{i,j} = ijMLI.FileNames;
    end
    MLIres.RootPath{i} = ijMLI.RootPath;
end
clear inp1 inp2
