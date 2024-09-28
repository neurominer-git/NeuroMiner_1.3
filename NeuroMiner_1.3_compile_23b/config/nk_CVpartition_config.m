% =========================================================================
% FORMAT res = nk_CVpartition_config(res)
% =========================================================================
% Setup repeated nested cross-validation structure 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 12/2023

function act = nk_CVpartition_config(defaultsfl, act)

global NM

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end
if ~exist('act','var') || isempty(act), act = 0; end

if ~defaultsfl

    CV2frm = 1;
    CV1ps = 'not defined'; CV1fs = CV1ps; CV2ps = CV1fs; CV2fs = CV2ps; 
    CV1pn = 10; CV1fn = 10; CV2pn = 10; CV2fn = 10; 
    Eq.enabled = false;
    Eq.bincount = 10;
    Eq.mincount = 10;
    Eq.maxcount = 10;
    Eq.addremoved2test = 1;
    Eq.Covar = NM.label;
    Eq.posnegrat = 1.5;
    Decomp = '';
    ConstrainedCV = 2;
    ConstrainedGroupIndex = [];
    CVOnlyGroupFlag = 2;
    CVOnlyGroupIndex = [];
    CVOnlyGroupLevel = 1;
    CV2ps = ['Define no. of CV2 permutations [ P2 = ' CV2ps ' ]|'];
    CV1ps = ['Define no. of CV1 permutations [ P1 = ' CV1ps ' ]|'];
    
    % Check whether cv structure already exists
    if isfield(NM,'cv')
        if iscell(NM.cv)
            nk_PrintLogo
            ncv = length(NM.cv);
            fprintf('\nMultiple CV structures detected:')
            fprintf('\n================================')
            for i=1:ncv
                % determine size of partitions
                fprintf('\n%g:\t CV2: [%g, %g], CV1: [%g, %g]', ...
                    i, size(NM.cv{i}.TrainInd,1), size(NM.cv{i}.TrainInd,2), ...
                    size(NM.cv{i}.cvin{1,1}.TrainInd,1), size(NM.cv{i}.cvin{1,1}.TrainInd,2))
            end
            actcv = nk_input('What to do',0,'m','Modify existing CV|Remove existing CV|Add new CV|<< Back',1:4);
            switch actcv
                case 1
                    selcv = nk_input(['Select CV for modification [1-' num2str(ncv) ']'],0,e);
                    tcv = NM.cv{selcv};
                case 2
                    selcv = nk_input(['Select CV for removal [1-' num2str(ncv) ']'],0,e);
                    tcv = NM.cv; NM = rmfield(NM,'cv'); cnt=1;
                    for i=1:ncv
                        if i==selcv, continue, end
                        NM.cv{cnt} = tcv{i};
                        cnt=cnt+1;
                    end
                    NM = nk_CVpartition_config(NM);
                case 3
                    selcv = ncv+1; tcv = [];
                case 4
                    return
            end

        else
            ncv = 1;
        end
    end
    fl = true;
    
    if strcmp(NM.modeflag, 'classification') && isfield(NM,'groupnames') && length(NM.groupnames) > 2
        Decomps = 'not defined';
        if isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND') && isfield(NM.TrainParam.RAND,'Decompose')
            switch NM.TrainParam.RAND.Decompose
                case 1
                    Decomps = 'One-vs-One';
                case 2
                    Decomps = 'One-vs-All';
                case 9
                    Decomps = 'Multi-group (no decomposition)';
            end
            Decompn = NM.TrainParam.RAND.Decompose;
        else
            Decompn = 1; fl = false;
        end
        Decomp = ['Define decomposition mode: ' Decomps '|'];
        MenuVec = [11,1:5];
    else
        NM.TrainParam.RAND.Decompose = 1; Decompn = 1;
        MenuVec = [11,1:4];
    end
    if isfield(NM.TrainParam.RAND,'Eq')
        Eq = NM.TrainParam.RAND.Eq;
    else
        NM.TrainParam.RAND.Eq = Eq;
    end
    if isfield(NM.TrainParam.RAND,'ConstrainedCV')
        ConstrainedCV = NM.TrainParam.RAND.ConstrainedCV;
        ConstrainedGroupIndex = NM.TrainParam.RAND.ConstrainedGroupIndex;
    else
        NM.TrainParam.RAND.ConstrainedCV = ConstrainedCV;
        NM.TrainParam.RAND.ConstrainedGroupIndex = ConstrainedGroupIndex;
    end

    if isfield(NM.TrainParam.RAND,'CVOnlyGroup')
        CVOnlyGroupFlag = NM.TrainParam.RAND.CVOnlyGroup.flag;
        CVOnlyGroupIndex = NM.TrainParam.RAND.CVOnlyGroup.index;
        CVOnlyGroupLevel = NM.TrainParam.RAND.CVOnlyGroup.level;
    else
        NM.TrainParam.RAND.CVOnlyGroup.flag = CVOnlyGroupFlag;
        NM.TrainParam.RAND.CVOnlyGroup.index = CVOnlyGroupIndex;
        NM.TrainParam.RAND.CVOnlyGroup.level = CVOnlyGroupLevel;
    end

    buildstr = ''; savestr = ''; loadstr = ''; MenuRem = []; CV2prx = '' ; CV1prx = ''; 
    EQstr = ''; EQcovstr = ''; EQminstr = ''; EQmaxstr = ''; EQbinstr = ''; EQshufflestr = ''; EQposnegstr=''; EQorigstr = ''; EQeqstr = '';
    constrainstr = ''; constraingrp = '';  yesno = {'yes','no'}; 
    cvgroupindexstr = ''; cvgrouplevelstr = ''; cvflagstr ='';

    %% Define menu options for cross-validation setup
    if isfield(NM,'TrainParam')

        if isfield(NM.TrainParam,'RAND')
            
            if isfield(NM.TrainParam.RAND,'lgoflp'), lgoflp = NM.TrainParam.RAND.lgoflp; end
            CV2reps = ''; MenuVec_extra = [];
            CV2STR_FRAME = '(Pooled) cross-validation';
            if isfield(NM.TrainParam.RAND,'CV2Frame')
               CV2frm = NM.TrainParam.RAND.CV2Frame;
               
               switch NM.TrainParam.RAND.CV2Frame
                    case 2
                       if NM.TrainParam.RAND.OuterPerm>1
                           CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Pooled ';
                       else
                           CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Pooled';
                       end
                       CV2reps = sprintf('Define number of Leave-Group-Out repetitions [ %g repetitions ]|',NM.TrainParam.RAND.OuterPerm);  MenuVec_extra = 12;
                    case 3
                       CV2STR_FRAME = 'Nested Leave-Group-Out';
                       NM.TrainParam.RAND.OuterPerm = 1;
                   case 4
                       CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Leave-Group-In';
                end
            end
            CV2Frame = ['Select cross-validation framework [ ' CV2STR_FRAME ' ]|'];

            if isfield(NM.TrainParam.RAND,'CV2Frame') && NM.TrainParam.RAND.CV2Frame ~= 1
                if ~isfield(NM.TrainParam.RAND,'CV2LCO')
                   CV2STR_LCO = 'not defined';
                   CV2fn = 1;
                else
                   CV2fn = numel(unique(NM.TrainParam.RAND.CV2LCO.ind));
                   CV2STR_LCO = sprintf('%g groups among %g cases',CV2fn, numel(NM.TrainParam.RAND.CV2LCO.ind));
                end
                MenuRem = 3;           
                CV2ps = ['Define Outer (CV2) Leave-Group-Out partitioning [ ' CV2STR_LCO ' ]|'];
                CV2fs = '';
            else
                if isfield(NM.TrainParam.RAND,'OuterPerm')
                    if isfield(NM.TrainParam.RAND,'OuterFold') && ...
                            (NM.TrainParam.RAND.OuterFold == -1 || ...
                            NM.TrainParam.RAND.OuterFold == numel(NM.label))
                        CV2ps = ''; CV2prx = ' [ LOO ]';
                        if strcmp(NM.modeflag,'classification'), NM.SVM.GridParam = 1; end
                        MenuRem = 2;
                    else
                        CV2pn = NM.TrainParam.RAND.OuterPerm; CV2ps = num2str(CV2pn);
                        CV2ps = ['Define no. of Outer (CV2) permutations [ P2 = ' CV2ps ' ]|'];
                    end
                else
                    fl = false;
                end
                if isfield(NM.TrainParam.RAND,'OuterFold')
                    CV2fn = NM.TrainParam.RAND.OuterFold; CV2fs = num2str(CV2fn);
                else
                    fl = false;
                end
                CV2fs = [ 'Define no. of Outer (CV2) folds [ K2 = ' CV2fs CV2prx ' ]|' ];
            end

            if isfield(NM.TrainParam.RAND,'InnerPerm')
                if isfield(NM.TrainParam.RAND,'CV2Frame') && NM.TrainParam.RAND.CV2Frame >= 3 
                    if ~isfield(NM.TrainParam.RAND,'CV1LCO')
                       CV1STR_LCO = 'not defined';
                       CV1fn = 1;
                    else
                       CV1fn = numel(unique(NM.TrainParam.RAND.CV1LCO.ind));
                       CV1STR_LCO = sprintf('%g groups among %g cases',CV1fn, numel(NM.TrainParam.RAND.CV1LCO.ind));
                    end
                    MenuRem = [3 5];    
                    switch NM.TrainParam.RAND.CV2Frame
                        case 3
                            CV1ps = ['Define Inner (CV1) Leave-Group-Out partitioning [ ' CV1STR_LCO ' ]|'];
                        case 4
                            CV1ps = ['Define Inner (CV1) Leave-Group-In partitioning [ ' CV1STR_LCO ' ]|'];
                    end
                    CV1fs = '';
                else
                    if isfield(NM.TrainParam.RAND,'InnerFold') && ... 
                            (NM.TrainParam.RAND.InnerFold == -1 || ...
                            NM.TrainParam.RAND.InnerFold == numel(NM.label) - floor(numel(NM.label) / CV2fn))
                        CV1ps = ''; CV1prx = ' [ LOO ]';
                        if strcmp(NM.modeflag,'classification'), NM.SVM.GridParam = 1; end
                        MenuRem = [MenuRem 4];
                    else
                        CV1pn = NM.TrainParam.RAND.InnerPerm; CV1ps = num2str(CV1pn);
                        CV1ps = ['Define no. of Inner (CV1) permutations [ P1 = ' CV1ps ' ]|'];
                    end
                    CV1fs = ['Define no. of Inner (CV1) folds [ K1 = ' num2str(NM.TrainParam.RAND.InnerFold) CV1prx ' ]|' ];
                end
            else
                fl = false;
            end
            if ~isfield(NM.TrainParam.RAND,'InnerFold') 
                fl = false;
            end
            if fl
                enabled = 2; enabledstr = 'no';
                addremoved2test = 2; addremoved2teststr = ''; shufflestr={'yes','no'};
                if isfield(NM.TrainParam,'RAND') && isfield(NM.TrainParam.RAND,'Eq') 
                    if isfield(NM.TrainParam.RAND.Eq,'enabled') && NM.TrainParam.RAND.Eq.enabled ==1, enabled = 1; enabledstr = 'yes'; end
                    if enabled == 1 && (isfield(NM.TrainParam.RAND.Eq,'addremoved2test') && NM.TrainParam.RAND.Eq.addremoved2test==1)
                        addremoved2test = 1; addremoved2teststr = ', shuffle to CV1 test data'; 
                    end
                end
                MenuVec = [MenuVec 13];
                switch NM.modeflag
                    case 'regression'
                        EQstr = ['Equalize label histogram at the CV1 cycle by undersampling [ ' enabledstr addremoved2teststr ' ]|']; 
                        if isfield(NM.TrainParam.RAND, 'Eq') && NM.TrainParam.RAND.Eq.enabled
                            EQcovstr = sprintf('Define target label(s) for uniform distribution [ size: %g,%g ]|', size(Eq.Covar));
                            EQminstr = sprintf('Define minimum # of observations at the lower end of label histogram [ %g ]|', Eq.mincount);
                            EQmaxstr = sprintf('Define minimum # of observations at the upper end of label histogram [ %g ]|', Eq.maxcount);
                            EQbinstr = sprintf('Define # of bins in label histogram [ %g ]|', Eq.bincount);
                            EQshufflestr = sprintf('Shuffle removed observations to CV1 test data [ %s ]|', shufflestr{addremoved2test});
                            EQorigstr = 'Show original histogram of the entire data|';
                            EQeqstr = 'Show (equalized) histogram of CV1 training partition [1,1] (build CV structure first!)|';
                            MenuVec = [MenuVec 14:20];
                        end
                    case 'classification'
                        EQstr = ['Equalize class sizes at the CV1 cycle by undersampling [ ' enabledstr addremoved2teststr ' ]|']; 
                        if isfield(NM.TrainParam.RAND,'Eq') && NM.TrainParam.RAND.Eq.enabled
                            EQposnegstr  = sprintf('Define positive/negative ratio after equalization (1=>same amount of positive/negative cases) [ %g ]|', Eq.posnegrat);
                            EQshufflestr = sprintf('Shuffle removed observations to CV1 test data [ %s ]|', shufflestr{addremoved2test});
                            MenuVec = [MenuVec 17 18 ];
                        end
                end
                if isfield(NM.TrainParam.RAND,'CV2Frame') && NM.TrainParam.RAND.CV2Frame == 1 && ~strcmp(CV2prx,' [ LOO ]') && ~strcmp(CV1prx,' [ LOO ]')
                    Consstr = yesno{ConstrainedCV};
                    constrainstr = ['Constrain Cross-Validation structure based on a group index variable [ ' Consstr ' ]|'];
                    MenuVec = [MenuVec 21];
                    if ConstrainedCV == 1
                        if ~isempty(ConstrainedGroupIndex)
                            constrindstr = sprintf('%g groups among %g cases', numel(unique(ConstrainedGroupIndex)), numel(ConstrainedGroupIndex));
                        else
                            constrindstr = 'undefined';
                        end
                        constraingrp = ['Define group index variable [ ' constrindstr ' ]|'];
                        MenuVec = [MenuVec 22];
                    end
                else
                    NM.TrainParam.RAND.ConstrainedCV=2;
                end
                
                CVflstr = yesno{CVOnlyGroupFlag};
                cvflagstr = ['Use part of the sample as exclusive test data [ ' CVflstr ']|']; MenuVec = [MenuVec 23];
                if CVOnlyGroupFlag == 1
                    if ~isempty(CVOnlyGroupIndex)
                        if islogical(CVOnlyGroupIndex)
                            cvgrpidxstr = sprintf('%g cases',sum(CVOnlyGroupIndex));
                        else
                            cvgrpidxstr = sprintf('%g cases',numel(CVOnlyGroupIndex));
                        end
                    else
                        cvgrpidxstr = 'undefined';
                    end
                    cvgroupindexstr = sprintf('Provide index vector (logical or numerical) test-only cases [ %s ]|',cvgrpidxstr); MenuVec = [MenuVec 24];
                    grpstr = {'CV2 test data', 'CV1 & CV2 test data'};
                    cvgrplvlstr = grpstr{CVOnlyGroupLevel};
                    cvgrouplevelstr = sprintf('TO BE IMPLEMENTED: Define cross-validation level where the exclusive test data should be used [ %s ]|', cvgrplvlstr); MenuVec = [MenuVec 25];
                end
                buildstr = 'Build Cross-Validation structure|';
                loadstr  = 'Load Cross-Validation structure|';
                MenuVec = [MenuVec 6 7];
                if isfield(NM,'cv')
                    savestr  = 'Save Cross-Validation structure|';
                    MenuVec  = [MenuVec 8];
                end
            end
        end
    end

    if ~act
        if ~isempty(CV2reps), MenuVec = [MenuVec(1:3) MenuVec_extra MenuVec(4:end) ]; end
        if ~isempty(MenuRem), MenuVec(MenuRem) = []; end
    
        nk_PrintLogo
    
        fprintf('\n****************************************************************************************')
        fprintf('\nSelect appropriate Kx (x=1/2) of folds for your dataset, where:')
        fprintf('\n\tKx < # of available subjects defines K-fold repeated, stratified cross-validation')
        fprintf('\n\tKx = -1 defines LOO cross-validation')
        fprintf('\nIf you choose Kx = -1 no CV permutations will be available.')
        fprintf('\n****************************************************************************************\n')
    
        act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> CROSS-VALIDATION SETTINGS',0,'mq', ...
                    [CV2Frame ...
                     CV2ps ...
                     CV2fs ...
                     CV2reps ...
                     CV1ps ...
                     CV1fs ...
                     Decomp ...
                     EQstr ...
                     EQcovstr ...
                     EQminstr ...
                     EQmaxstr ...
                     EQbinstr ...
                     EQposnegstr ...
                     EQshufflestr ...
                     EQorigstr ...
                     EQeqstr ...
                     constrainstr ...
                     constraingrp ...
                     cvflagstr ...
                     cvgroupindexstr ...
                     cvgrouplevelstr ...
                     buildstr loadstr savestr], MenuVec,1);
    end
    switch act
       
        case 1
             if isfield(NM.TrainParam.RAND,'CV2Frame') 
                 switch NM.TrainParam.RAND.CV2Frame
                     case 1
                        NM.TrainParam.RAND.OuterPerm = ...
                            nk_input('Number of permutations for Outer (CV2) cross-validation',0,'w1',CV2pn);
                     case {2,3,4}
                         NM.TrainParam.RAND.CV2LCO.ind = ...
                             nk_input('Define index vector for Outer (CV2) case-to-group assignment',0,'i',[],[numel(NM.cases),1]);
                         if NM.TrainParam.RAND.CV2Frame~=2
                             NM.TrainParam.RAND.OuterPerm = 1;
                         end
                         NM.TrainParam.RAND.OuterFold = numel(unique(NM.TrainParam.RAND.CV2LCO.ind));
                 end
             else
                 NM.TrainParam.RAND.OuterPerm = ...
                            nk_input('Number of permutations for Outer (CV2) cross-validation ',0,'w1',CV2pn);
             end
        case 2
             OuterFold = nk_input('Number of folds for Outer (CV2) cross-validation ',0,'i',CV2fn);
             if OuterFold <= 1 && OuterFold ~= -1
                 errordlg('Specify at least 2 outer (CV2) folds for k-fold cross-validation and 1 for training-only mode!');
             else
                 NM.TrainParam.RAND.OuterFold = OuterFold;
             end
        case 3
             if isfield(NM.TrainParam.RAND,'CV2Frame') 
                 switch NM.TrainParam.RAND.CV2Frame
                    case {1,2}
                        NM.TrainParam.RAND.InnerPerm = ...
                            nk_input('Number of permutations for Inner (CV1) cross-validation',0,'w1',CV1pn);
                    case {3,4}    
                          NM.TrainParam.RAND.CV1LCO.ind = ...
                                 nk_input('Define index vector for Inner (CV1) case-to-group assignment',0,'i',[],[numel(NM.cases),1]);
                             NM.TrainParam.RAND.InnerPerm = 1;
                             NM.TrainParam.RAND.InnerFold = numel(unique(NM.TrainParam.RAND.CV1LCO.ind));
                 end
             else
                  NM.TrainParam.RAND.InnerPerm = ...
                            nk_input('Number of permutations for Inner (CV1) cross-validation',0,'w1',CV1pn);
             end
             
        case 4
             InnerFold = nk_input('Number of folds for Inner cross-validation (CV1)',0,'i',CV1fn);
             if InnerFold < 1 && InnerFold ~= -1
                 errordlg('Specify at least 2 inner (CV1) folds for k-fold cross-validation and 1 for training-only mode!');
             else
                 NM.TrainParam.RAND.InnerFold = InnerFold;
             end

        case 5
             if NM.TrainParam.RAND.InnerFold == -1 || NM.TrainParam.RAND.OuterFold == -1
                 NM.TrainParam.RAND.Decompose = ...
                     nk_input('Multi-group decomposition method',0,'m', ...
                     'One-vs-All',2);
             else
                 NM.TrainParam.RAND.Decompose = ...
                     nk_input('Multi-group decomposition method',0,'m', ...
                     'One-vs-One|One-vs-All|Multi-Class [not fully implemented yet]',[1,2,9],Decompn);
             end
         
        case 6
             groups = []; appendfl=false; oldcv=zeros(0,size(NM.label,2));

             % Check whether to overwrite or append to current CV structure
             if isfield(NM,'cv') && ~isempty(NM.cv)
                 if NM.TrainParam.RAND.OuterFold == size(NM.cv(1).TrainInd,2) && ...
                    NM.TrainParam.RAND.InnerFold == size(NM.cv(1).cvin{1,1}.TrainInd,2)
                        appendfl = nk_input('Append to existing CV structure',0,'m', ...
                            ['No|' ...
                            'Append to Outer CV cycle|' ...
                            'Append to Inner CV cycles|' ...
                            'Rebuild CV structure without new permutations'],0:3,1);
                        if appendfl, oldcv = NM.cv; end
                 end
             end
             if isfield(NM,'groupnames')
                 groupnames = NM.groupnames;
             else
                 groupnames = [];
             end
             if strcmp(NM.modeflag,'classification')
                 for i=1:size(NM.label,2)
                     if ~isempty(oldcv), ioldcv = oldcv(i); else, ioldcv=[]; end
                     if size(NM.label,2)>1, grpn = groupnames{i}; else, grpn = groupnames; end
                     cv(i) =  nk_MakeCrossFolds(NM.label(:,i), NM.TrainParam.RAND, NM.modeflag, groups, grpn, ioldcv, appendfl);
                 end
             else
                 cv = nk_MakeCrossFolds(NM.label, NM.TrainParam.RAND, NM.modeflag, groups, groupnames, oldcv, appendfl);
             end
             if ~isempty(cv)
                 for i=1:numel(cv)
                     switch appendfl
                         case 1
                            NM.cv(i).TrainInd = [NM.cv(i).TrainInd; cv(i).TrainInd];
                            NM.cv(i).TestInd 	= [NM.cv(i).TestInd; cv(i).TestInd];
                            NM.cv(i).cvin	= [NM.cv(i).cvin; cv(i).cvin];
                            if strcmp(NM.modeflag,'classification')
                                NM.cv(i).class = [NM.cv(i).class; cv(i).class];
                                NM.cv(i).classnew = [NM.cv(i).classnew; cv(i).classnew];
                            end
                         otherwise
                             if ~isfield(NM,'cv') || isempty(NM.cv)
                                NM.cv = cv(i);
                             else
                                if ~isfield(NM,'class') && i==1
                                    NM = rmfield(NM,'cv');
                                end
                                NM.cv(i) = cv(i);
                             end
                     end
                 end
             end
        case 7
             loadcv = nk_FileSelector(1,'matrix','Select cross-validation structure', 'CVstruct_.*\.mat');
             if ~isempty(loadcv)
                load(loadcv,'cv');
                NM.cv = cv;
             end
             % Note CV (24.05.2023): 
             % to do: update RAND: 
             % RAND: CV2Frame, OuterFold, OuterPerm, InnerFold, InnerPerm,
             % Decompose, Eq (what information can we get from 'cv', what
             % else do we need?) 
        case 8
             savecv = nk_input('Filename for cross-validation structure (will be saved in current workind directory)',0,'s');
             if ~isempty(savecv)
                 cv =  NM.cv;
                 save(['CVstruct_' savecv '.mat'],'cv');
             end
        case 11
             NM.TrainParam.RAND.CV2Frame = nk_input('Select cross-validation framework',0,'m', ... 
                 'Pooled|Outer Leave-Group-Out/Inner pooled|Nested Leave-Group-Out|Outer Leave-Group-Out/Inner Leave-Group-In',1:4,CV2frm);
             switch NM.TrainParam.RAND.CV2Frame
                 case 1
                    if isfield(NM.TrainParam.RAND,'CV2LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV2LCO'); end
                    if isfield(NM.TrainParam.RAND,'CV1LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV1LCO'); end
                 case 2
                    if isfield(NM.TrainParam.RAND,'CV1LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV1LCO'); end
             end
         case 12
             NM.TrainParam.RAND.OuterPerm = nk_input('Number of repetitions for Outer (CV2) cross-validation',0,'w1',CV2pn);
         case 13
             if ~isfield(NM.TrainParam.RAND,'Eq')
                 NM.TrainParam.RAND.Eq.enabled = false;
             end
             NM.TrainParam.RAND.Eq.enabled = ~NM.TrainParam.RAND.Eq.enabled ;
         case 14
             if isfield(NM,'covars') && ~isempty(NM.covars)
                eqtarget = nk_input('Perform histogram equalization using the target label or some other covariate', ...
                                    0,'m','Target label|Covariate',[0,1],1);
                if eqtarget
                     covind = nk_SelectCovariateIndex(NM);
                     NM.TrainParam.RAND.Eq.Covar = NM.covars(:,covind);
                     NM.TrainParam.RAND.Eq.CovarName = NM.covnames{covind};
                else
                    if isfield(NM.TrainParam,'MULTILABEL')
                         NM.TrainParam.RAND.Eq.Covar = nm_nanmean(NM.label(:,NM.TrainParam.MULTILABEL.sel),2);
                     else
                         NM.TrainParam.RAND.Eq.Covar = nm_nanmean(NM.label,2);
                     end
                     NM.TrainParam.RAND.Eq.CovarName  = 'NM.label';
                end
             end
        case 15
             NM.TrainParam.RAND.Eq.mincount = ...
                             nk_input('Minimum # of observation at the lower end of label histogram',0,'e', Eq.mincount);
        case 16    
             NM.TrainParam.RAND.Eq.maxcount = ...
                 nk_input('Minimum # of observation at the upper end of label histogram',0,'e', Eq.maxcount);
        case 17
            switch NM.modeflag
                 case 'classification'
                     NM.TrainParam.RAND.Eq.posnegrat = nk_input('Target ratio bigger / smaller class (1 => classes have sample size after balancing)',0,'e', Eq.posnegrat);
                 case 'regression'
                      NM.TrainParam.RAND.Eq.bincount = nk_input('# Bins for histogram analysis',0,'i',Eq.bincount);
            end
        case 18
             NM.TrainParam.RAND.Eq.addremoved2test = ~NM.TrainParam.RAND.Eq.addremoved2test;
        case 19
            figure;histogram(NM.TrainParam.RAND.Eq.Covar, NM.TrainParam.RAND.Eq.bincount); 
            ylabel('# of observations in bins'); 
            xlabel(NM.TrainParam.RAND.Eq.CovarName); 
            title(['Histogram analysis of ' Eq.CovarName]);
        case 20 
            figure; histogram(NM.TrainParam.RAND.Eq.Covar(NM.cv.TrainInd{1,1}(NM.cv.cvin{1,1}.TrainInd{1,1})), NM.TrainParam.RAND.Eq.bincount);
            ylabel('# of observations in bins'); 
            xlabel(NM.TrainParam.RAND.Eq.CovarName); 
            title(['Analysis of equalized histogram of ' Eq.CovarName]);
        case 21
            if NM.TrainParam.RAND.ConstrainedCV == 1, NM.TrainParam.RAND.ConstrainedCV = 2; else, NM.TrainParam.RAND.ConstrainedCV = 1 ; end
        case 22
            NM.TrainParam.RAND.ConstrainedGroupIndex = ...
                                 nk_input('Define index vector for constrained, stratified CV',0,'i',[],[numel(NM.cases),1]);
        case 23
            if NM.TrainParam.RAND.CVOnlyGroup.flag == 1, NM.TrainParam.RAND.CVOnlyGroup.flag = 2; else, NM.TrainParam.RAND.CVOnlyGroup.flag = 1; end
        case 24
            NM.TrainParam.RAND.CVOnlyGroup.index = ...
                                 nk_input('Define index vector for exclusive test data use of cases',0,'i',[],[numel(NM.cases),1]);
        case 25
            if NM.TrainParam.RAND.CVOnlyGroup.level == 1, NM.TrainParam.RAND.CVOnlyGroup.level = 2; else, NM.TrainParam.RAND.CVOnlyGroup.level = 1; end

    end
else
    NM.TrainParam.RAND.CV2Frame     = 1;
    NM.TrainParam.RAND.OuterFold    = 10;
    NM.TrainParam.RAND.OuterPerm    = 1;
    NM.TrainParam.RAND.InnerFold    = 10;
    NM.TrainParam.RAND.InnerPerm    = 1;
end

