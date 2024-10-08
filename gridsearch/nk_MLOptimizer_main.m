function analysis = nk_MLOptimizer_main(inp, datid, PreML)

global SVM GRD MODEFL W2AVAIL VERBOSE BATCH

if ~isfield(inp,'batchflag') || isempty(inp.batchflag), inp.batchflag = false; end
if ~isfield(inp,'probflag') || isempty(inp.probflag), inp.probflag = false; end
if ~isfield(inp,'simFlag') || isempty(inp.simFlag), inp.simFlag = false; end
BATCH = inp.batchflag;
[ ~, strout] = GetMLType(SVM);

W2AVAIL = nk_CheckW2Avail(SVM, MODEFL);

if VERBOSE, fprintf('\n************ %s-TRAINING ************', strout); end

% ********************* Setup of input parameters *************************
% Create input parameter structure for nk_MLOptimizer
if numel(inp.X) > 1
    inparams.Y = inp.X(inp.curmodal).Y; 
    if isfield(inp.X(inp.curmodal),'Yw'), inparams.Yw = inp.X(inp.curmodal).Yw; end
else
    inparams.Y = inp.X.Y; 
    if isfield(inp.X,'Yw'), inparams.Yw = inp.X.Yw; end
end

inparams.id =               datid;
inparams.stranalysis =      inp.stranalysis;
inparams.nclass  =          inp.nclass;
inparams.multiflag =        inp.multiflag;
inparams.probflag =         inp.probflag;
inparams.batchflag =        inp.batchflag;
inparams.GridAct =          inp.GridAct;
inparams.ovrwrtGD =         inp.ovrwrtGD;
inparams.updGD =            inp.updGD;
inparams.ovrwrtRes =        inp.ovrwrtRes;
inparams.updRes =           inp.updRes;
inparams.label =            inp.labels;
inparams.l =                size(inp.labels,1);
inparams.simFlag =          inp.simFlag;
if strcmp(MODEFL,'classification')
    inparams.ngroups =      numel(unique(inp.labels(~isnan(inp.labels))));           
else
    inparams.ngroups =      1;
end
inparams.covars =           inp.covars;
inparams.covars_oocv =      [];
inparams.varstr =           inp.varstr;
inparams.varind =           inp.tF;
inparams.P =                inp;
inparams.analyses =         inp.analyses;
inparams.stacking =         inp.stacking;
inparams.rootdir =          inp.rootdir;
inparams.maindir =          inp.maindir;
inparams.gdmat =            inp.gdmat;
inparams.preprocmat =       inp.preprocmat;
if isfield(inp,'ParamOptimizationMode')
    inparams.ParamOptimizationMode = inp.ParamOptimizationMode;
end
if isfield(inp,'time2event')
    inparams.time2event =   inp.time2event;
end
if isfield(inp, 'C')
    inparams.C =            inp.C;
end
inparams.curanal =          inp.curanal;
clear inp

% Prepare Parameter containers
inparams.Params = cell(inparams.nclass,1);
inparams.Params_desc = cell(inparams.nclass,1);
inparams.nMLparams = zeros(inparams.nclass,1);
inparams.nPreMLparams = zeros(inparams.nclass,1);
inparams.ModalityVec = cell(inparams.nclass,1);

   
% Define ML optimization params
for i=1:inparams.nclass        
    
    %% Define hyperparameters (C, nu, epsilon) of SVM provided by LIBSVM
    % (this piece of code originates from a development state of NM where
    % only LIBSVM was provided as algorithm (it will be deleted once one
    % common interface for all algorithms provided by NM is implemented)
    switch SVM.prog 
        
        case 'MEXELM'
            inparams.Params{i}{end+1} = GRD.Cparams;
            inparams.Params_desc{i}{end+1} = 'Slack';
      
        case {'LIBSVM','LIBLIN'}
            % Check classifier type (L1/L2 or Nu / One-class)
            ctype = nk_GetLIBSVMClassType(SVM);
            % Check if regression parameters are needed 
            rtype = nk_GetLIBSVMRegrType(SVM);
            switch rtype
                case 1
                    inparams.Params{i}{1} = GRD.Epsparams;
                    inparams.Params_desc{i}{1} = 'Epsilon';
                case 2
                    inparams.Params{i}{1} = GRD.Nuparams;
                    inparams.Params_desc{i}{1} = 'Nu';
            end
            % This is the slack / nu-SVC parameter of the SVM
            switch ctype 
                case 1
                    inparams.Params{i}{end+1} = GRD.Nuparams;
                    inparams.Params_desc{i}{end+1} = 'Nu (SVC)';
                otherwise
                    inparams.Params{i}{end+1} = GRD.Cparams;
                    inparams.Params_desc{i}{end+1} = 'Slack';
            end
            
        otherwise
            if VERBOSE, fprintf('\n%s #%g: no slack parameters needed.',strout,i); end
    end

    %% Define kernel parameter settings
    switch SVM.kernel.kernstr
        case {' -t 0', 'lin','linear','lin_kernel','none','lin_elm'}
            fprintf('\n%s #%g: no kernel parameters needed.',strout,i)     
        case {' -t 1', 'poly', 'polynomial', 'Polynomial', 'polyN', 'hpolyN'}
            inparams.Params{i}{end+1} = GRD.Gparams;
            inparams.Params_desc{i}{end+1} = 'Kernel';
            inparams.Params{i}{end+1} = GRD.PolyCoefparams;
            inparams.Params_desc{i}{end+1} = 'Poly coef';
            inparams.Params{i}{end+1} = GRD.PolyDegrparams;
            inparams.Params_desc{i}{end+1} = 'Poly Degree';
        case ' -t 2'
            inparams.Params{i}{end+1} = GRD.Gparams;
            inparams.Params_desc{i}{end+1} = 'Kernel';
        case ' -t 3'
            inparams.Params{i}{end+1} = GRD.Gparams;
            inparams.Params_desc{i}{end+1} = 'Kernel';
            inparams.Params{i}{end+1} = GRD.PolyCoefparams;
            inparams.Params_desc{i}{end+1} = 'Sigmoid coef';   
        case ' -t 4'
            if SVM.kernel.kerndef < 8
                inparams.Params{i}{end+1} = GRD.WLiters;
                inparams.Params_desc{i}{end+1} = 'Kernel';
            else
                if isfield(SVM.kernel,'customfunc_nargin')
                    if SVM.kernel.customfunc_nargin >= 1
                        for n = 1:SVM.kernel.customfunc_nargin
                            argName = sprintf('customkernel_arg%d', n); 
                            %GRD.(argName) = 0;
                            inparams.Params{i}{end+1} = GRD.(argName);
                            inparams.Params_desc{i}{end+1} = sprintf('Kernel function argument %d', n);
                        end
                    end
                end
            end
        otherwise
            switch SVM.prog
                case 'LIBSVM'
                    switch SVM.LIBSVM.LIBSVMver
                        case 3
                            switch SVM.kernel.kernstr
                                case {' -t 4', ' -t 5'}
                                    inparams.Params{i}{end+1} = GRD.PolyCoefparams;
                                    inparams.Params_desc{i}{end+1} = 'Kernel';
                                case {' -t 6', ' -t 7'}
                                     inparams.Params{i}{end+1} = GRD.Gparams;
                                     inparams.Params_desc{i}{end+1} = 'Kernel';
                                otherwise
                                    error('Precomputed kernels are not supported in this version of NM');
                            end
                        otherwise
                            error('Precomputed kernels are not supported in this version of NM');
                    end
                otherwise
                    % This is the gamma exponent parameter of the RBF kernel
                    inparams.Params{i}{end+1} = GRD.Gparams;
                    inparams.Params_desc{i}{end+1} = 'Kernel';
            end
    end
    
    %% other algorithms have to be treated as special cases (will be revised in the future to have one parameter interface for all algorithms)
    switch SVM.prog
        case {'matLRN','GLMNET','GRDBST','ROBSVM','ELASVM'}
            if strcmp(SVM.prog,'matLRN'), F='matLearn'; else, F=SVM.prog; end
            if ~isfield(GRD.(F),'Params') || isempty(GRD.(F).Params)
                 if VERBOSE, fprintf('\n%s #%g: No parameters have to be optimized.',strout,i); end
            else
                for j=1:numel(GRD.(F).Params)
                    inparams.Params{i}{end+1} = GRD.(F).Params(j).range;
                    inparams.Params_desc{i}{end+1} = GRD.(F).Params(j).name;
                end
            end

        case 'MLPERC'
            for j=1:numel(GRD.(SVM.prog).Params)
                %Use the index for 'MLPERC_hidden_layer_sizes'
                if strcmp(GRD.(SVM.prog).Params(j).name, 'MLPERC_hidden_layer_sizes')
                    inparams.Params{i}{end+1} = 1:size(GRD.(SVM.prog).Params(j).range, 1);
                else
                    inparams.Params{i}{end+1} = GRD.(SVM.prog).Params(j).range;
                end
                
                inparams.Params_desc{i}{end+1} = GRD.(SVM.prog).Params(j).name;
            end

        case 'MEXELM'
            inparams.Params{i}{end+1} = GRD.Neuronparams;
            inparams.Params_desc{i}{end+1} = 'Hidden neurons';
        case 'kNNMEX'
            inparams.Params{i}{end+1} = GRD.Kparams;
            inparams.Params_desc{i}{end+1} = 'Nearest neighbors';
        case 'BLOREG'
            inparams.Params{i}{end+1} = GRD.Tolparams;
            inparams.Params_desc{i}{end+1} = 'Tolerance';
        case {'LIBSVM','LIBLIN','CCSSVM'}
            if SVM.(SVM.prog).Weighting
                if isfield(GRD,'Weightparams')
                    inparams.Params{i}{end+1} = GRD.Weightparams;
                    inparams.Params_desc{i}{end+1} = 'Weight Factor';
                else
                    inparams.Params{i}{end+1} = 1;
                    inparams.Params_desc{i}{end+1} = 'Weight Factor';
                end
                if isfield(GRD,'CCLambdaparams')
                    inparams.Params{i}{end+1} = GRD.CCLambdaparams;
                    inparams.Params_desc{i}{end+1} = 'CC-Lambda';
                end
            end
        case 'DECTRE'
            inparams.Params{i}{1} = GRD.Leafparams;
            inparams.Params_desc{i}{1} = 'Leafness';
        case 'RNDFOR'
            inparams.Params{i}{1} = GRD.Treeparams;
            inparams.Params_desc{i}{1} = 'Decision Trees';
            inparams.Params{i}{2} = GRD.NumDparams;
            inparams.Params_desc{i}{2} = 'Num Feats';
            inparams.Params{i}{3} = GRD.RFMaxDparams;
            inparams.Params_desc{i}{3} = 'Criterion';
            inparams.Params{i}{4} = GRD.RFCritdefsparams;
            inparams.Params_desc{i}{4} = 'Max Depth';
            inparams.Params{i}{5} = GRD.RFMinSSplitparams;
            inparams.Params_desc{i}{5} = 'Min Samples Split';
            inparams.Params{i}{6} = GRD.RFMinSLeafparams;
            inparams.Params_desc{i}{6} = 'Min Samples Leaf';
            inparams.Params{i}{7} = GRD.RFMinWFLeafparams;
            inparams.Params_desc{i}{7} = 'Min Weight Fract Leaf';
            inparams.Params{i}{8} = GRD.RFMaxLeafNparams;
            inparams.Params_desc{i}{8} = 'Max Leaf Nodes';
            inparams.Params{i}{9} = GRD.RFMinImpDecrparams;
            inparams.Params_desc{i}{9} = 'Min Impurity Decrease';
            inparams.Params{i}{10} = GRD.RFBootstrapparams;
            inparams.Params_desc{i}{10} = 'Bootstrap';
            inparams.Params{i}{11} = GRD.RFOobScoreparams;
            inparams.Params_desc{i}{11} = 'Out-Of-Bag Score';
            inparams.Params{i}{12} = GRD.RFCcpAlphaparams;
            inparams.Params_desc{i}{12} = 'CCP Alpha';
            inparams.Params{i}{13} = GRD.RFMaxSampparams;
            inparams.Params_desc{i}{13} = 'Max Sample';
            if strcmp(MODEFL, 'classification')
                inparams.Params{i}{14} = GRD.RFClassWeightparams;
                inparams.Params_desc{i}{14} = 'Class Weigths';
            end

        case 'SEQOPT'
            inparams.Params{i}{1} = 1:size(SVM.SEQOPT.C,1);
            inparams.Params_desc{i}{1} = 'Sequence';
            inparams.Params{i}{2} = GRD.CutOffparams;
            inparams.Params_desc{i}{2} = '%Step';
            inparams.Params{i}{3} = GRD.LimsLparams;
            inparams.Params_desc{i}{3} = '%L-Population';
            inparams.Params{i}{4} = GRD.LimsUparams;
            inparams.Params_desc{i}{4} = '%U-Population';
        case 'IMRELF'
            inparams.Params{i}{1} = GRD.Cparams;
            inparams.Params_desc{i}{1} = 'Lambda';
            inparams.Params{i}{2} = GRD.Gparams;
            inparams.Params_desc{i}{2} = 'Sigma';
        case 'WBLCOX'
            inparams.Params{i}{1} = GRD.CoxCutoffparams;
            inparams.Params_desc{i}{1} = 'P-Cox-cutoff';
    end

    % Define number of ML algo parameter to be optimized in
    % the nk_GridSearch process
    inparams.nMLparams(i) = numel(inparams.Params{i});

    % Add further parameters to parameter vector if they had been used
    % during data-preprocessing
    if ~isempty(PreML) 
        inparams.nPreMLparams(i) = numel(PreML.Params_desc);
        for p = 1:inparams.nPreMLparams(i)
            inparams.Params_desc{i}{end+1} = PreML.Params_desc{p};
            inparams.Params{i}{end+1} = PreML.Params{p};
        end
        inparams.ModalityVec{i} = PreML.ModalityVec;
    end
end

analysis = nk_MLOptimizer(inparams, inparams.stranalysis, inparams.id, inparams.GridAct, inparams.batchflag);

end
