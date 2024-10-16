function [PREPROC, act, stepind] = nk_Preproc_config(PREPROC, varind, parentstr, stepind, enind)
% =========================================================================
% FORMAT [PREPROC, act, stepind] = nk_Preproc_config(PREPROC, varind, ...
%                                                parentstr, stepind, enind)
% =========================================================================
% Configuration of Data Preprocessing Pipelines:
% This is a menu-driven NM wizard for setting up pre-processing pipelines  
% to be applied to data modalities prior to running the proper machine 
% learning algorithms. 
%
% Pre-CV steps: Data operations are performed on a case-by-case basis and
% therefore only need to be run once - before the data are split according
% to NM.cv
%
% In-CV steps: Data operations are dependent on group-level processing and
% results and therefore have to be fully wrapped into the rdCV framework
%
% Post-processing steps: Data operation are performed after the completion
% of the pre-processing pipeline in the respective CV1 partition (currently
% only label imputation in the training samples).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 09/2024

% Defaults:
% ---------
global NM EXPERT DEV

% Set Step Index to 1 if none is specified
if ~exist('stepind','var') || isempty(stepind), stepind = 1; end

if isempty(PREPROC)
    fprintf('No preprocessing pipeline defined for modality %g', varind)
    copyflag = nk_input(['Do you want to use another pipeline as template for modality #' num2str(varind) ' definitions?'],0,'yes|no',[1,0],1);
    if copyflag
        varind_copy = nk_input('Specify modality index',0,'e');
        if varind_copy <= numel(NM.TrainParam.PREPROC) && ~isempty(NM.TrainParam.PREPROC{varind_copy})
            PREPROC = NM.TrainParam.PREPROC{varind_copy};
        end
    else
        nan_in_label = false; if sum(isnan(NM.label(:)))>0, nan_in_label=true; end
        nan_in_pred=false; if sum(isnan(NM.Y{varind}(:)))>0, nan_in_pred=true; end
        PREPROC = nk_DefPREPROC_config(NM.modeflag, nan_in_pred, nan_in_label);
    end
end

if isfield(NM.TrainParam,'LABEL') && NM.TrainParam.LABEL.flag == 1
    modeflag = NM.TrainParam.LABEL.newmode;
else
    modeflag = NM.modeflag;
end

if ~exist('enind','var'), enind = []; end

imaging_flag = false;
if NM.datadescriptor{varind}.type == 1
    imaging_flag = true;
end

if ~isstruct(enind)
    fl = true;
    if (EXPERT || DEV) && imaging_flag
        EF = struct('correctnuis', fl, ...
            'reducedim',     fl, ...
            'extdim',        fl, ...
            'standardize',   fl, ...
            'scale',         fl, ...
            'normalize',     fl, ...
            'unitnormalize', fl, ...
            'binning',       fl, ...
            'impute',        fl, ...
            'labelimpute',   fl, ...
            'elimzero',      fl, ...
            'remmeandiff',   fl, ...
            'rankfeat',      fl, ...
            'remvarcomp',    fl, ...
            'devmap',        fl, ...
            'JuSpace',       fl,...
            'ROImeans',      fl);
    elseif (EXPERT || DEV) && ~imaging_flag
        EF = struct('correctnuis', fl, ...
            'reducedim',     fl, ...
            'extdim',        fl, ...
            'standardize',   fl, ...
            'scale',         fl, ...
            'normalize',     fl, ...
            'unitnormalize', fl, ...
            'binning',       fl, ...
            'impute',        fl, ...
            'labelimpute',   fl, ...
            'elimzero',      fl, ...
            'remmeandiff',   fl, ...
            'rankfeat',      fl, ...
            'remvarcomp',    fl, ...
            'devmap',        fl);

        if DEV
            EF.graphComputation = fl;
            EF.graphSparsity = fl;
            EF.graphMetrics = fl;
            EF.customPreproc = fl;
        end     
    elseif imaging_flag 
        EF = struct('correctnuis', fl, ...
            'reducedim',     fl, ...
            'extdim',        fl, ...
            'standardize',   fl, ...
            'scale',         fl, ...
            'normalize',     fl, ...
            'unitnormalize', fl, ...
            'binning',       fl, ...
            'impute',        fl, ...
            'labelimpute',   fl, ...
            'elimzero',      fl, ...
            'remmeandiff',   fl, ...
            'rankfeat',      fl, ...
            'remvarcomp',    fl, ...
            'devmap',        fl,...
            'JuSpace',       fl,...
            'ROImeans',      fl);
    else
        EF = struct('correctnuis', fl, ...
            'reducedim',     fl, ...
            'extdim',        fl, ...
            'standardize',   fl, ...
            'scale',         fl, ...
            'normalize',     fl, ...
            'unitnormalize', fl, ...
            'binning',       fl, ...
            'impute',        fl, ...
            'labelimpute',   fl, ...
            'elimzero',      fl, ...
            'remmeandiff',   fl, ...
            'rankfeat',      fl, ...
            'remvarcomp',    fl, ...
            'devmap',        fl);
    end

else
    EF = enind;
end
actstr=[]; actmnu=[]; prestr=[];

% Get description of parameters
d = nk_GetParamDescription2(NM, PREPROC, 'PreProc');
nk_PrintLogo

% Check group processing mode possibilities
if isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
    label_temp = NM.TrainParam.LABEL.newlabel; 
else
    label_temp = NM.label; 
end

if numel(unique(label_temp)) > 2 && strcmp(modeflag,'classification') 
    if isfield(NM.TrainParam,'RAND') && ...
                isfield(NM.TrainParam.RAND,'Decompose') && ...
                    NM.TrainParam.RAND.Decompose == 2
        fprintf('\nOne-vs-All mode => multigroup processing activated');
        PREPROC.BINMOD = 0;
    else
        fprintf('\n%s',d.PREPROC.groupmode);
        cmdstr = 'Define group processing mode in multi-class setting'; cmdmnu = 1;
        [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu); 
    end
elseif isempty(PREPROC) || ~isfield(PREPROC,'BINMOD')
    PREPROC = config_binmod(NM, PREPROC);
end
 
% Check for target scaling / exponential transformation option
if strcmp(modeflag,'regression')
    if isfield(PREPROC,'LABELMOD') && isfield(PREPROC.LABELMOD,'TARGETSCALE') && ( PREPROC.LABELMOD.TARGETSCALE || isfield(PREPROC.LABELMOD,'POLYNOM') )
        cmdstr = 'Modify / Disable';
        if ~strcmp(d.PREPROC.targetscaling,'NA') 
            prestr = sprintf('\n* %s', d.PREPROC.targetscaling); 
        end
    else
        cmdstr = 'Enable';
    end
    cmdstr = [cmdstr ' label transformation']; cmdmnu = 99;
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
end

% Check for availability of image filtering options
spatfltflag = false; imganalflag = false;
if NM.datadescriptor{varind}.type == 1 && NM.TrainParam.STACKING.flag==2
    imganalflag = true;
    if isfield(PREPROC,'SPATIAL') && PREPROC.SPATIAL.cubetype>1
        cmdstr = 'Modify / Disable ';
        if ~strcmp(d.PREPROC.spatialfiltering,'NA') 
            prestr = sprintf('%s\n* %s', prestr, d.PREPROC.spatialfiltering); 
        end
        spatfltflag = true;
    else
        cmdstr = 'Enable';
    end
    cmdstr = [cmdstr ' spatial operations using Spatial OP Wizard ']; cmdmnu = 98;
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
end

if ~isempty(prestr)
    fprintf('\n'); fprintf('NON-CV PREPROCESSING STEPS \n');
    fprintf('========================== ')
    fprintf('%s ',prestr)
end

slnan = sum(isnan(label_temp));
if slnan
    fprintf('\n');
    cmdstr = 'Define parameters for label propagation to unlabeled training data'; cmdmnu = 100;
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    if ~strcmp(d.PREPROC.labelimpute,'NA')
        fprintf('Missing label settings:\n* %s', d.PREPROC.labelimpute); 
    else
        fprintf('Missing labels found! Please specifiy label imputation parameters')
    end
end

%% Setup main sequence menu
if isfield(PREPROC,'ACTPARAM') && ~isempty(PREPROC.ACTPARAM)
    
    lact = length(PREPROC.ACTPARAM);
    
    %% Display current preprocessing sequence
    fprintf('\n\n')
    if NM.TrainParam.STACKING.flag==2
        fprintf('CV-PREPROCESSING PIPELINE \n')
        fprintf('========================= ')
    else
        fprintf('CV-PREPROCESSING PIPELINE [STACKING MODE] \n')
        fprintf('========================================= ')
    end
    dimredflag = false;
    for i=1:numel(d.PREPROC.preprocact)
        stepstr = sprintf('Step %g: %s',i,d.PREPROC.preprocact{i});
        switch PREPROC.ACTPARAM{i}.cmd
            case 'extfeat'
                if i==stepind+1
                    fprintf('\n\t'); fprintf('\\__ %s ',stepstr);
                else
                    fprintf('\n\t    %s ',stepstr);
                end
            case 'extdim'
                if i==stepind
                    fprintf('\n\t'); fprintf('\\__ %s ',stepstr);
                else
                    fprintf('\n\t    %s',stepstr);
                end
            case 'reducedim'
                dimredflag=true;
                if i==stepind || (i == stepind-1 && strcmp(PREPROC.ACTPARAM{i+1}.cmd,'extdim'))
                    fprintf('\n'); fprintf('>> %s ',stepstr);                 
                else
                    fprintf('\n   %s', stepstr);
                end
            case {'scale','standardize','normalize','unitnormalize'}
                if i==stepind
                    if imganalflag && ~spatfltflag && ~dimredflag
                        fprintf('\n'); fprintf('>> %s [!!! CV1/2 test data offset errors expected without spatial smoothing !!!]',stepstr); 
                    else
                        fprintf('\n'); fprintf('>> %s ',stepstr); 
                    end
                else
                    if imganalflag && ~spatfltflag && ~dimredflag
                        fprintf('\n   %s [!!! CV1/2 test data offset errors expected without spatial smoothing !!!]',stepstr); 
                    else
                        fprintf('\n   %s',stepstr); 
                    end
                end
            case 'remvarcomp'
                if ~isfield(PREPROC.ACTPARAM{i}.REMVARCOMP,'G')
                    fprintf('\n'); fprintf('   %s ',stepstr);
                else
                    if i==stepind 
                        fprintf('\n'); fprintf('>> %s ',stepstr); 
                    else
                        fprintf('\n   %s ', stepstr);
                    end
                end
            otherwise
                if i==stepind 
                    fprintf('\n'); fprintf('>> %s ',stepstr); 
                else
                    fprintf('\n   %s ', stepstr);
                end
        end
    end
    
    cmdstr = 'Add preprocessing step|Remove preprocessing step';
    cmdmnu = [2 3];
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    
    if lact > 1 
        
        if stepind == lact
            cmdstr = ['Insert preprocessing step' ...
                      '|Replace current preprocessing step' ...
                      '|Modify current preprocessing step' ...
                      '|<< Previous step'];
            cmdmnu = [4 5 6 7];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
             
        elseif stepind == 1
            cmdstr = ['Insert preprocessing step' ...
                      '|Replace current preprocessing step' ...
                      '|Modify current preprocessing step' ...
                      '|Next step >>'];
            cmdmnu = [4 5 6 8];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        else
            cmdstr = ['Insert preprocessing step' ...
                      '|Replace current preprocessing step' ...
                      '|Modify current preprocessing step' ...
                      '|Next step >>' ...
                      '|<< Previous step'];
            cmdmnu = [4 5 6 8 7];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        end
        cmdstr = 'Change order of preprocessing steps';
        cmdmnu = 9;
        [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        if lact > 3
            cmdstr = 'Go to step ...'; cmdmnu = 10;
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        end
    else
        cmdstr = ['Insert preprocessing step' ...
                  '|Replace current preprocessing step' ...
                  '|Modify current preprocessing step' ];
        cmdmnu = [4 5 6];
        [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    end
    
    titlestr = sprintf('Preprocessing pipeline generator: %g / %g steps.', stepind, lact);

else
    cmdstr = 'Add preprocessing step';
    cmdmnu = 2;
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    titlestr = sprintf('Preprocessing pipeline generator: no steps defined.');
end

mestr = titlestr; navistr = [parentstr ' >>> ' mestr]; fprintf('\n\nYou are here: %s >>> ',parentstr); 
act = nk_input(titlestr, 0,'mq', actstr, actmnu);

switch act
    
    case 98
        PREPROC = config_spatialfilter(PREPROC, navistr, NM.brainmask{varind});
        
    case 99
        PREPROC = config_targetscaling(PREPROC, navistr);
        
    case 100
        if ~isfield(PREPROC,'LABELMOD')
            PREPROC.LABELMOD=[];
        end
        PREPROC.LABELMOD = config_labelimpute(PREPROC.LABELMOD, navistr);
        
    case 1 % Configure group processing mode
        PREPROC = config_binmod(NM, label_temp, PREPROC);
    
    case 2 % Add Preprocessing step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, [], 0, EF, navistr, modeflag);
        if isfield(PREPROC,'ACTPARAM')
            if strcmp(PREPROC.ACTPARAM{end}.cmd, 'extfeat')
                stepind = numel(PREPROC.ACTPARAM)-1;
            elseif strcmp(PREPROC.ACTPARAM{end}.cmd, 'extdim') 
                stepind = numel(PREPROC.ACTPARAM)-1;
            else
                stepind = numel(PREPROC.ACTPARAM);
            end
        end
    case 3 % Remove preprocessing step
        if lact > 1
            if strcmp(PREPROC.ACTPARAM{stepind}.cmd, 'rankfeat')
                % Delete also extfeat
                remind = [stepind stepind+1];
            elseif strcmp(PREPROC.ACTPARAM{stepind}.cmd, 'reducedim') && numel(PREPROC.ACTPARAM) > stepind && strcmp(PREPROC.ACTPARAM{stepind+1}.cmd, 'extdim')
                remind = [stepind stepind+1];
            else
                remind = stepind;
            end
            PREPROC.ACTPARAM(remind) = [];
        else
            PREPROC = rmfield(PREPROC,'ACTPARAM');
        end
        if stepind > 1
            if any(strcmp(PREPROC.ACTPARAM{stepind-1}.cmd,{'extfeat','extdim'})) && stepind - 2 > 0
                stepind = stepind - 2;
            else
                stepind = stepind - 1;
            end
        end
    case 4 % Insert Preprocessing step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 2, EF, navistr, modeflag);
    
    case 5 % Replace current step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 1, EF, navistr, modeflag);
        
    case 6 % Modify current step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 0, EF, navistr, modeflag);
    
    case 7 % Go to previous step
        if strcmp(PREPROC.ACTPARAM{stepind-1}.cmd,'extfeat') && stepind - 2 > 0
            stepind = stepind - 2;
        else
            stepind = stepind - 1;
        end
    
    case 8 % Go to next step
        if strcmp(PREPROC.ACTPARAM{stepind}.cmd,'rankfeat') && ...
                stepind + 2 <= numel(PREPROC.ACTPARAM)
            stepind = stepind + 2;
%         elseif strcmp(PREPROC.ACTPARAM{stepind}.cmd,'reducedim') && ...
%                 (stepind + 2 <= numel(PREPROC.ACTPARAM) && strcmp(PREPROC.ACTPARAM{stepind+1}.cmd,'extdim'))  
%             stepind = stepind + 2;
        else
            stepind = stepind + 1;
        end
        
    case 9 % Change order of preprocessing steps
        neworder = nk_input('Define new order of preprocessing steps',0,'e',[],[1, lact]);
        fl = true(1,2); newcmdorder = cell(numel(neworder),1); 
        for i = 1:numel(neworder)
            newcmdorder{i} = PREPROC.ACTPARAM{neworder(i)}.cmd; 
        end
        for i = 1:numel(neworder)-1
            if strcmp(newcmdorder{i},'rankfeat') && ~strcmp(newcmdorder{i+1},'extfeat') 
                fl(1) = false;  
            end
            if any(strcmp({'reducedim','remvarcomp'}, newcmdorder{i})) && ~strcmp(newcmdorder{i+1},'extdim')
                for j=1:numel(neworder)
                    if strcmp(PREPROC.ACTPARAM{j}.cmd,'extdim'), fl(2) = false; break; end 
                end
            end
        end
        
        if sum(fl) == numel(fl)
            PREPROC.ACTPARAM = PREPROC.ACTPARAM(neworder);
        else
            mess = 'Reordering cannot be performed!'; 
            if ~fl(1)
                mess = char(mess,'Weighting-based feature generation has to be executed right after the computation of the feature weighting / ranking!');
            end
            if ~fl(2)
                mess = char(mess,'Extraction of dimensions has to be executed right after the computation of the dimensionality reduction model!');
            end
            mess = cellstr(mess);
            mess = sprintf('%s\n',mess{:});
            errordlg(mess);
        end
        stepind = 1;
        
    case 10
        tstepind = nk_input('Go to preprocessing step',0,'w1',stepind);
        if tstepind > numel(PREPROC.ACTPARAM)
            tstepind = numel(PREPROC.ACTPARAM); 
        elseif tstepind < 1
            tstepind = 1;
        end
        if strcmp(PREPROC.ACTPARAM{tstepind}.cmd,'extfeat') 
            stepind = tstepind -1;
        elseif any(strcmp({'reducedim','remvarcomp'}, PREPROC.ACTPARAM{tstepind}.cmd)) && (tstepind+1 == numel(PREPROC.ACTPARAM) && strcmp(PREPROC.ACTPARAM{tstepind+1}.cmd, 'extdim'))
            stepind = tstepind -1;
        else
            stepind = tstepind;
        end

end

if ~isfield(PREPROC,'BINMOD'), PREPROC.BINMOD = 1; end

end

% -------------------------------------------------------------------------
function [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, ...
                                                replflag, EF, navistr, modeflag)

if ~isempty(PREPROC) && isfield(PREPROC,'ACTPARAM'), ...
    lact = length(PREPROC.ACTPARAM); 
else
    lact = 0; 
end
if ~exist('stepind','var') || isempty(stepind)
    modflag = 0; stepind = lact + 1;
else
    if replflag == 2, modflag = 0; else, modflag = 1; end
end

M = NM.TrainParam.FUSION.M;
fusefl = NM.TrainParam.FUSION.flag;
                    
if ~exist('replflag','var') || isempty(replflag), replflag = 0; end

actstr = []; actmnu = []; CURACT2 = [];

if modflag && ~replflag
    
    CURACT = PREPROC.ACTPARAM{stepind};
    cmdstr = CURACT.cmd;
    switch cmdstr
        case 'impute'
            cmd = 2;
        case 'correctnuis'
            cmd = 3;
        case 'scale'
            cmd = 4;
        case 'normalize'
            cmd = 5;
        case 'standardize'
            cmd = 6;
        case {'discretize','symbolize'}
            cmd = 7;
        case 'reducedim'
            cmd = 8;
        case 'labelimpute'
            cmd = 9;
        case 'elimzero'
            cmd = 10;
        case 'rankfeat'
            cmd = 11;
            CURACT2 = PREPROC.ACTPARAM{stepind+1};
        case 'remmeandiff'
            cmd = 12;
        case 'unitnormalize'
            cmd = 13;
        case 'remvarcomp'
            cmd = 14;
        case 'extdim'
            cmd = 15;
        case 'devmap'
            cmd = 16;
        case 'graphSparsity'
            cmd = 17;
        case 'graphMetrics'
            cmd = 18;
        case 'graphComputation'
            cmd = 19;
        case 'JuSpace'
            cmd = 20;
        case 'ROImeans'
            cmd = 21;
        case 'customPreproc'
            cmd = 22;
    end
   
else    
    fld = fieldnames(EF);
    for z = 1:numel(fld)
        if EF.(fld{z})
            switch fld{z}
                 case 'impute'
                    switch fusefl 
                        case 1
                            tY = cell2mat(NM.Y(M));
                        otherwise
                            tY = NM.Y{varind};
                    end
                    nans = isnan(tY); clear tY; 
                    if sum(nans(:)) || NM.TrainParam.STACKING.flag==1
                        if NM.TrainParam.STACKING.flag==2
                            snan_rows = sum(nans,2); snan_mat = sum(snan_rows);
                            cmdstr = sprintf('Impute missing values (%g NaNs in %1.1f%% of cases) ', ...
                                snan_mat, sum(snan_rows>0)*100/size(nans,1));                           cmdmnu = 2;
                        else
                            cmdstr = sprintf('Impute missing input-layer predictions');                 cmdmnu = 2;
                        end
                    else
                        continue
                    end
                case 'correctnuis'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Regress out nuisance covariates'; cmdmnu = 3;
                    else
                        continue;
                    end
                case 'scale'
                    cmdstr = 'Scale data';                                                          cmdmnu = 4;
                case 'normalize'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Normalize to group mean';                                         cmdmnu = 5;
                    else
                        continue;
                    end
                case 'standardize'
                    cmdstr = 'Standardize data';                                                    cmdmnu = 6;
                case 'binning'
                    cmdstr = 'Apply binning method to data';                                        cmdmnu = 7;
                case 'reducedim'
                    cmdstr = 'Apply dimensionality reduction method to data';                       cmdmnu = 8;
                case 'labelimpute'
                    if isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
                        slnan = sum(isnan(NM.TrainParam.LABEL.newlabel));
                    else
                        slnan = sum(isnan(NM.label));
                    end
                    if slnan
                        cmdstr = 'Propagate labels to unlabeled training data';                     cmdmnu = 9;
                    else
                        continue
                    end
                case 'elimzero'
                    cmdstr = 'Prune non-informative columns from data matrix';                      cmdmnu = 10;
                case 'rankfeat'
                    cmdstr = 'Rank / Weight features';                                              cmdmnu = 11;
                case 'remmeandiff'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Remove group-level differences using offset correction';          cmdmnu = 12;
                    else
                        continue;
                    end
                case 'unitnormalize'
                    cmdstr = 'Normalize data to unit vector';                                       cmdmnu = 13;
                case 'remvarcomp'
                    cmdstr = 'Extract variance components from data';                               cmdmnu = 14;
                case 'extdim'
                    if stepind > 1 && any(strcmp({'reducedim','remvarcomp'}, PREPROC.ACTPARAM{stepind-1}.cmd))
                        cmdstr = 'Extract subspaces from reduced data projections';                 cmdmnu = 15;
                    end
                case 'devmap'
                    cmdstr = 'Measure deviation from normative data';                               cmdmnu = 16;
                case 'graphSparsity'
                    cmdstr = 'Apply sparsity threshold to connectivity matrix';                     cmdmnu = 17;
                case 'graphMetrics'
                    cmdstr = 'Compute network metrics from connectivity matrices';                  cmdmnu = 18;
                case 'graphComputation'
                    cmdstr = 'Compute individual networks from input data';                         cmdmnu = 19;
                case 'JuSpace'
                    cmdstr = 'Correlation with neurotransmitter systems (PET oder SPECT maps; JuSpace Toolbox)';                cmdmnu = 20; 
                case 'ROImeans'
                    cmdstr = 'Compute ROI mean values';                                             cmdmnu = 21; 
                case 'customPreproc'
                    cmdstr = 'Add a custom preproc function (from .m-file)';                        cmdmnu = 22;

            end
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
            cmdstr =[]; cmdmnu=[];
        end
    end   
    nk_PrintLogo
    cmd = nk_input(['Select preprocessing type [Step ' ...
        num2str(stepind) ' / ' num2str(lact) ']'],0,'mq',actstr,actmnu);
    if ~cmd, return; end 
    CURACT = [];
   
end

switch cmd
    case 2
        CURACT = config_impute(NM, varind, CURACT, navistr);
    case 3
        CURACT = config_covars(NM, varind, CURACT, navistr);
    case 4
        CURACT = config_scaling(CURACT, navistr);
    case 5
        CURACT = config_groupnorm(NM, varind, CURACT, navistr);
    case 6
        CURACT = config_standard(CURACT, navistr);
    case 7
        CURACT = config_binning(CURACT);
    case 8
        CURACT = config_dimred(CURACT, navistr);
    case 9
        CURACT = config_labelimpute(CURACT, navistr);
    case 10
        CURACT = config_elimzero(CURACT, navistr);
    case 11
        CURACT = config_rankfeat( NM, varind, CURACT, navistr );
        CURACT2 = config_waction( NM, varind, CURACT2, navistr );
    case 12
        CURACT = config_remmeandiff(NM, CURACT);
    case 13
        CURACT = config_unitnorm( CURACT, navistr );
    case 14
        CURACT = config_remvarcomp( NM, varind, CURACT, navistr );
    case 15
        switch  PREPROC.ACTPARAM{stepind-1}.cmd
            case 'reducedim'
                CURACT = config_extdim( CURACT, PREPROC.ACTPARAM{stepind-1}.DR , navistr );
            case 'remvarcomp'
                CURACT = config_extdim( CURACT, PREPROC.ACTPARAM{stepind-1}.REMVARCOMP, navistr );
        end
    case 16
        CURACT = config_devmap( NM, CURACT, navistr );
    case 17
        CURACT = config_graphSparsity(CURACT, navistr);
    case 18
        CURACT = config_graphMetrics(CURACT, navistr);
    case 19
        CURACT = config_graphConstruction(NM, varind, CURACT, navistr);
    case 20
        CURACT = config_JuSpace(NM, varind, CURACT, navistr);
    case 21
        CURACT = config_ROImeans(NM, varind, CURACT, navistr);
    case 22
        CURACT = config_customPreproc(CURACT, navistr);
    
end

switch replflag
    case 2 % Insert
        tACTPARAM = PREPROC.ACTPARAM(stepind:end);
        stepper = 1; if ~isempty(CURACT2), stepper = 2; end
        PREPROC.ACTPARAM(stepind+stepper:end+stepper) = tACTPARAM;
end
if ~isempty(CURACT),PREPROC.ACTPARAM{stepind} = CURACT; end
if ~isempty(CURACT2),PREPROC.ACTPARAM{stepind+1} = CURACT2; end

% Check if template preprocessing has to be activated 
if isfield(PREPROC,'TEMPLPROC'), PREPROC = rmfield(PREPROC,'TEMPLPROC'); end
for i=1:lact
    if isfield(PREPROC.ACTPARAM{i},'TEMPLPROC') && PREPROC.ACTPARAM{i}.TEMPLPROC
        PREPROC.TEMPLPROC = PREPROC.ACTPARAM{i}.TEMPLPROC; 
        break
    end
end
if isfield(PREPROC,'USEALL'), PREPROC = rmfield(PREPROC,'USEALL'); end
for i=1:lact
    if isfield(PREPROC.ACTPARAM{i},'USEALL') && PREPROC.ACTPARAM{i}.USEALL
        PREPROC.USEALL = PREPROC.ACTPARAM{i}.USEALL; 
        break
    end
end

end
% -------------------------------------------------------------------------
%%%% Target scaling (in regression) %%%%
function CURACT = config_targetscaling(CURACT, navistr)

if ~isfield(CURACT,'LABELMOD'), LABELMOD = []; else, LABELMOD = CURACT.LABELMOD; end
act = 1; while act > 0, [LABELMOD, act] = nk_LabelMod_config(LABELMOD, navistr ); end
CURACT.LABELMOD = LABELMOD;

end
% -------------------------------------------------------------------------
%%%% GROUP PROCESSING MODE (if multi-class available) %%%%
function CURACT = config_binmod(NM, label, CURACT)

% Default parameter
if isfield(CURACT,'BINMOD') && ~isempty(CURACT.BINMOD)
    switch CURACT.BINMOD, case 1, tBINMOD = 1; case 0, tBINMOD = 2; case 2, tBINMOD = 3; end
else
    if max(label(:,1))>2 && strcmp(modeflag,'classification') 
        tBINMOD = 0; % Multi-group
    else
        tBINMOD = 1; % Binary mode
    end
end

if max(label(:,1))<=2 || ~strcmp(NM.modeflag,'classification')
    CURACT.BINMOD = 1;
else
    CURACT.BINMOD = nk_input('Group processing mode',0, 'm', 'binary|multi-group',[1,0], tBINMOD);
end

end

% -------------------------------------------------------------------------
%%%% IMPUTATION CONFIGURATION %%%%
function CURACT = config_impute(NM, varind, CURACT, navistr)

if ~isfield(CURACT,'IMPUTE'), CURACT.IMPUTE = []; end
act = 1; while act > 0, [CURACT.IMPUTE, act] = nk_Impute_config(NM, CURACT.IMPUTE, varind, navistr); end
CURACT.cmd = 'impute';

end

% -------------------------------------------------------------------------
%%%% LABEL IMPUTATION CONFIGURATION %%%%
function CURACT = config_labelimpute(CURACT, navistr)

if ~isfield(CURACT,'LABELIMPUTE'), CURACT.LABELIMPUTE = []; end
act = 1; while act > 0, [CURACT.LABELIMPUTE, act] = nk_LabelImpute_config(CURACT.LABELIMPUTE, navistr); end
CURACT.cmd = 'labelimpute';

end

% -------------------------------------------------------------------------
%%%% NUISANCE COVARIATES %%%%
function CURACT = config_covars(NM, varind, CURACT, navistr)

if ~isfield(CURACT,'PX'), CURACT.PX = []; end

act = 1; while act > 0, [CURACT, act] = nk_PartialCorrelations_config(NM, CURACT, varind, navistr); end

CURACT.cmd = 'correctnuis';

end

% -------------------------------------------------------------------------
%%%% CORRECT FOR MULTIPLCATIVE EFFECTS %%%%
function CURACT = config_groupnorm(NM, varind, CURACT, navistr)

if isfield(NM,'covars') && ~isempty(NM.covars)
    act = 1; while act > 0, [CURACT, act] = nk_GroupNorm_config(NM, varind, CURACT, navistr); end
    CURACT.cmd = 'normalize';
end

end
% -------------------------------------------------------------------------
%%%% CORRECT FOR GROUP-WISE OFFSETS %%%%
function CURACT = config_remmeandiff(NM, CURACT)

tSIND = 'NA';
tDIND = 'NA';

if isfield(NM,'covars') && ~isempty(NM.covars)
    
    if ~isempty(CURACT) && isfield(CURACT,'IND') && ( ~isempty(CURACT.IND) || ~strcmp(CURACT.IND,'NA'))
        tSIND = CURACT.sIND;
    else
        tSIND = 1;
    end

    tSIND = nk_input('Select group index covariate for the computation of global and group-specific means (zeros indicate cases to be skipped)',0,'i', tSIND);
    
    if ~isempty(CURACT) && isfield(CURACT,'dIND') && ( ~isempty(CURACT.dIND) || ~strcmp(CURACT.dIND,'NA'))
        tDIND = CURACT.dIND;
    else
        tDIND = 1;
    end

    tDIND = nk_input('Select group index covariate to apply global and group-specific means corrections to (zeros indicate cases to be skipped)',0,'i', tDIND);
    
end

CURACT.sIND = tSIND;
CURACT.dIND = tDIND;
CURACT.cmd = 'remmeandiff';

end

% -------------------------------------------------------------------------
%%%% SCALING %%%%
function CURACT = config_scaling( CURACT, navistr )
if ~isfield(CURACT,'SCALE'), CURACT.SCALE = []; end
act = 1; while act > 0, [CURACT.SCALE, act] = nk_Scale_config(CURACT.SCALE, navistr); end
CURACT.cmd = 'scale';

end

% -------------------------------------------------------------------------
%%%% SYMBOLIZATION / DISCRETIZATION %%%%
function CURACT = config_binning(CURACT)

act = 1; while act, [CURACT, act] = nk_DiscSym_config(CURACT, 'MAIN >> NM workspace configuration >> Preprocessing pipeline generator'); end

end

% -------------------------------------------------------------------------
%%%% STANDARDIZATION %%%%
function CURACT = config_standard(CURACT, navistr)

act = 1; while act > 0, [CURACT, act] = nk_Standardize_config( CURACT, navistr); end
CURACT.cmd = 'standardize';

end

% -------------------------------------------------------------------------
%%%% UNIT NORMALIZATION %%%%
function CURACT = config_unitnorm(CURACT, navistr)

if ~isfield(CURACT,'UNITNORM'), CURACT.UNITNORM = []; end
act = 1; while act >0, [ CURACT, act ] = nk_Unitnorm_config(CURACT, navistr); end
CURACT.cmd = 'unitnormalize';

end

% -------------------------------------------------------------------------
%%%% ELIMINATE Zero var | NaN | Inf ATTRIBUTES %%%%
function CURACT = config_elimzero(CURACT, navistr)

if ~isfield(CURACT,'PRUNE'), CURACT.PRUNE=[]; end
if ~isfield(CURACT,'PX'), CURACT.PX = []; end
act = 1; while act >0, [ CURACT.PRUNE, CURACT.PX, act ] = nk_Prune_config(CURACT.PRUNE, CURACT.PX, navistr); end
CURACT.cmd = 'elimzero';

end

% -------------------------------------------------------------------------
%%%% DIM. REDUCTION %%%%
function CURACT = config_dimred(CURACT, navistr)

if isfield(CURACT,'DR'); DR = CURACT.DR; else, DR = []; end
if isfield(CURACT,'PX'), PX = CURACT.PX ; else, PX = []; end
if isfield(CURACT, 'TEMPLPROC'), TEMPLPROC = CURACT.TEMPLPROC; else, TEMPLPROC = []; end
if isfield(CURACT, 'USEALL'), USEALL = CURACT.USEALL; else, USEALL = []; end
act = 1; while act >0, [ DR, PX, TEMPLPROC, USEALL, act ] = nk_DimRed_main_config(DR, PX, TEMPLPROC, USEALL, navistr); end
if ~isempty(DR), CURACT.DR = DR; CURACT.PX = PX; CURACT.TEMPLPROC = TEMPLPROC; CURACT.USEALL = USEALL; CURACT.cmd = 'reducedim'; end

end
% -------------------------------------------------------------------------
%%%% Extraction of dimensionalities %%%%
function CURACT = config_extdim(CURACT, DR, navistr)

if isfield(CURACT,'PX') 
    PercMode                = CURACT.EXTDIM.PercMode;
    RedMode                 = CURACT.EXTDIM.RedMode;
    dims                    = nk_ReturnParam('dimensions', CURACT.PX.Px.Params_desc, CURACT.PX.opt);
else
    if ~isfield(DR,'PercMode')
        PercMode                = DR.dimmode;
        RedMode                 = 'PCA';
    else
        PercMode                = DR.PercMode;
        RedMode                 = DR.RedMode;
    end
    
    CURACT.PX               = nk_AddParam(DR.dims, 'dimensions', 1, []); dims = DR.dims;
end
act = 1; while act > 0, [dims, PercMode, act] = nk_ExtDim_config(RedMode, PercMode, dims, 0, navistr); end
CURACT.EXTDIM.dims = dims; CURACT.EXTDIM.PercMode = PercMode; CURACT.EXTDIM.RedMode = RedMode;
CURACT.PX = nk_AddParam(CURACT.EXTDIM.dims, 'dimensions', 1, CURACT.PX);
CURACT.cmd = 'extdim';

end

% -------------------------------------------------------------------------
%%%% SPATIAL FILTERING %%%%
function CURACT = config_spatialfilter( CURACT , navistr, brainmask)

if isfield(CURACT,'SPATIAL') 
    SPATIAL = CURACT.SPATIAL;
    if isfield(CURACT.SPATIAL,'PX'), PX = CURACT.SPATIAL.PX; end
else
    SPATIAL = []; PX = [];
end
[SPATIAL, PX ] = nk_Spatial_config(SPATIAL, PX, [], navistr, brainmask);
if ~isempty(SPATIAL), SPATIAL.PX = PX; CURACT.SPATIAL = SPATIAL; end

end

% -------------------------------------------------------------------------
%%%% Feature ranking %%%%
function CURACT = config_rankfeat(NM, varind, CURACT, navistr)

if isfield(CURACT,'RANK'), RANK = CURACT.RANK; else, RANK = []; end
if isfield(CURACT,'PX'), PX = CURACT.PX ; else, PX = []; end
[CURACT.RANK, CURACT.PX] = nk_Rank_config(RANK, PX, NM, varind, [], navistr);
CURACT.cmd = 'rankfeat';

end

% -------------------------------------------------------------------------
%%%% Feature extraction based on feature ranking / weighting %%%%
function CURACT = config_waction( NM, varind, CURACT, navistr )

if isfield(CURACT,'W_ACT'), W_ACT = CURACT.W_ACT; else, W_ACT = []; end
datadesc = NM.datadescriptor{varind}; brainmask = [];
if datadesc.type, brainmask = NM.brainmask{varind}; end
if isfield(CURACT,'PX'), PX = CURACT.PX ; else, PX = []; end
[ CURACT.W_ACT, CURACT.PX ] = nk_WAction_config( W_ACT, PX, datadesc, brainmask, [], navistr ); 
CURACT.cmd = 'extfeat';

end

% -------------------------------------------------------------------------
%%%% Removal of variance components %%%%
function CURACT = config_remvarcomp( NM, varind, CURACT, navistr )

if isfield(CURACT,'REMVARCOMP'), REMVARCOMP = CURACT.REMVARCOMP; else, REMVARCOMP = []; end
if isfield(CURACT,'PX'), PX = CURACT.PX ; else, PX = []; end
[ CURACT.REMVARCOMP, CURACT.PX ] = nk_Remvarcomp_config( NM, varind, REMVARCOMP, PX, navistr );
CURACT.cmd = 'remvarcomp';

end

% -------------------------------------------------------------------------
%%%% Mapping of deviation from normative data %%%%
function CURACT = config_devmap( NM, CURACT, navistr )

if isfield(CURACT,'DEVMAP'), DEVMAP = CURACT.DEVMAP; else, DEVMAP = []; end
if isfield(CURACT,'PX'), PX = CURACT.PX ; else, PX = []; end
[ CURACT.DEVMAP, CURACT.PX ] = nk_Devmap_config( DEVMAP, PX, NM, [], navistr );
CURACT.cmd = 'devmap';

end

function CURACT = config_graphSparsity(CURACT, navistr)

if ~isfield(CURACT,'GRAPHSPARSITY'), CURACT.GRAPHSPARSITY=[]; end
if ~isfield(CURACT,'PX'), CURACT.PX = []; end
act = 1; while act >0, [ CURACT.GRAPHSPARSITY, CURACT.PX, act ] = cv_graphSparsity_config(CURACT.GRAPHSPARSITY, CURACT.PX, navistr); end
CURACT.cmd = 'graphSparsity';

end

function CURACT = config_graphMetrics(CURACT, navistr)

if ~isfield(CURACT,'GRAPHMETRICS'), CURACT.GRAPHMETRICS=[]; end
if ~isfield(CURACT,'PX'), CURACT.PX = []; end
act = 1; while act >0, [ CURACT.GRAPHMETRICS, CURACT.PX, act ] = cv_graphMetrics_config(CURACT.GRAPHMETRICS, CURACT.PX, navistr); end
CURACT.cmd = 'graphMetrics';

end

function CURACT = config_graphConstruction(NM, varind, CURACT, navistr)

if ~isfield(CURACT,'GRAPHCONSTRUCTION'), CURACT.GRAPHCONSTRUCTION=[]; end
if ~isfield(CURACT,'PX'), CURACT.PX = []; end
datadesc = NM.datadescriptor{varind}; brainmask = [];
if datadesc.type, brainmask = NM.brainmask{varind}; end
act = 1; while act >0, [ CURACT.GRAPHCONSTRUCTION, CURACT.PX, act ] = cv_graphConstruction_config(CURACT.GRAPHCONSTRUCTION, CURACT.PX, brainmask, navistr); end
CURACT.cmd = 'graphComputation';

end

% -------------------------------------------------------------------------
function CURACT = config_JuSpace(NM, varind, CURACT, navistr)

if ~isfield(CURACT,'JUSPACE'), CURACT.JUSPACE=[]; end
%if ~isfield(CURACT,'PX'), CURACT.PX = []; end
datadesc = NM.datadescriptor{varind}; brainmask = [];
if datadesc.type, brainmask = NM.brainmask{varind}; end
act = 1; while act >0, [ CURACT.JUSPACE, act ] = JuSpace_config(CURACT.JUSPACE, brainmask, navistr); end
CURACT.cmd = 'JuSpace';
		
end

% -------------------------------------------------------------------------
function CURACT = config_ROImeans(NM, varind, CURACT, navistr)

if ~isfield(CURACT,'ROIMEANS'), CURACT.ROIMEANS=[]; end
%if ~isfield(CURACT,'PX'), CURACT.PX = []; end
datadesc = NM.datadescriptor{varind}; brainmask = [];
if datadesc.type, brainmask = NM.brainmask{varind}; end
act = 1; while act >0, [ CURACT.ROIMEANS, act ] = cv_ROImeans_config(CURACT.ROIMEANS, brainmask, navistr); end
CURACT.cmd = 'ROImeans';

end

% -------------------------------------------------------------------------
function CURACT = config_customPreproc(CURACT, navistr)

if ~isfield(CURACT,'CUSTOMPREPROC'), CURACT.CUSTOMPREPROC=[]; end
if ~isfield(CURACT,'PX'), CURACT.PX = []; end
act = 1; while act >0, [ CURACT.CUSTOMPREPROC, CURACT.PX, act ] = cv_customPreproc_config(CURACT.CUSTOMPREPROC, CURACT.PX, navistr); end
CURACT.cmd = 'customPreproc';

end


% -------------------------------------------------------------------------
function [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu)

if isempty(actstr),
    actstr = cmdstr;
    actmnu = cmdmnu;
else
    actstr = [ actstr '|' cmdstr ];
    actmnu = [ actmnu cmdmnu ];
end

end
