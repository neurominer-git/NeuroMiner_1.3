function [ STANDARD, act ] = nk_Standardize_config(STANDARD, parentstr, defaultsfl)
global NM

methodarr = {'standardization using median','standardization using mean','mean-centering','l1-median centering','qn-standardization','sn-standardization'};
sIND = []; dIND = []; WINSOPT = []; IQRFUN = 1; ZeroOut = 1; CALIBUSE = 2; METHOD = methodarr{1};
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    %% Retrieve previous settings
    if exist('STANDARD','var') && ~isempty(STANDARD)
        if isfield(STANDARD,'PX'), PX = STANDARD.PX; end
        if isfield(STANDARD,'METHOD'), METHOD = STANDARD.METHOD; end
        if isfield(STANDARD,'CALIBUSE'), CALIBUSE = STANDARD.CALIBUSE; end
        if isfield(STANDARD,'WINSOPT'), WINSOPT = STANDARD.WINSOPT; end
        if isfield(STANDARD,'IQRFUN'), IQRFUN = STANDARD.IQRFUN; end
        if isfield(STANDARD,'zerooutflag'), ZeroOut = STANDARD.zerooutflag; end
        if isfield(STANDARD,'sIND') && ~isempty(STANDARD.sIND), sIND = STANDARD.sIND; end 
        if isfield(STANDARD,'dIND') && ~isempty(STANDARD.dIND), dIND = STANDARD.dIND; end
    end
    
    %% Set defaults and menu strings
    if isfield(STANDARD,'sIND') && ~isempty(STANDARD.sIND)
        if CALIBUSE == 1
            sINDSTR = 'yes, using calibration data'; 
        else
            sINDSTR = sprintf('yes, using training data [ %s ] ', strjoin(NM.covnames(STANDARD.sIND)));
        end
        sINDDEF = 1;
    else
        sINDSTR = 'no'; sINDDEF = 2;
    end
    if isfield(STANDARD,'dIND') && ~isempty(STANDARD.dIND)
        dINDSTR = sprintf('yes, [ %s ]', strjoin(NM.covnames(STANDARD.dIND))); dINDDEF = 1;
    else
        dINDSTR = 'no'; dINDDEF = 2;
    end
    if ~isempty(WINSOPT) 
        WINSOPTSTR = sprintf('yes, threshold: +/- %s SD',nk_ConcatParamstr(WINSOPT, true)); WINSOPTDEF = 1;
    else
        WINSOPTSTR = 'no'; WINSOPTDEF = 2;
    end
    if IQRFUN == 1
        IQRFUNSTR = 'yes'; 
    else
        IQRFUNSTR = 'no'; 
    end
    if ZeroOut == 1
        ZEROOUTSTR = 'yes'; 
    else
        ZEROOUTSTR = 'no';
    end
    
    menuact = [ 1 2 3];
    WINSOPTENTRY = []; IQRENTRY = [];
    if any(strcmp(METHOD,{'standardization using mean', 'standardization using median'}))
        WINSOPTENTRY = ['Winsorize data (clamping of outliers) [ ' WINSOPTSTR ' ]|' ];
        menuact = 1:4;
        if strcmp(METHOD,'standardization using median')
            IQRENTRY = ['Use Inter-Quartile Range (IQR) instead of standard deviation [ ' IQRFUNSTR ' ]|' ];
            menuact = [ menuact 5 ];
        else
            IQRENTRY = [];
        end
    end
    menuact = [menuact 6];
    
    menustr = ['Select standardization method [ ' METHOD ' ]|', ...
               'Compute standardization model uing a subgroup of cases [ ' sINDSTR ' ]|', ...
               'Apply standardization model to a subgroup of cases [ ' dINDSTR ' ]|', ...
               WINSOPTENTRY, ...
               IQRENTRY, ...
               'Zero-out completely non-finite features [ ' ZEROOUTSTR ' ]'];
    
    [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);
    
    %% Display menu
    nk_PrintLogo
    mestr = 'Standardization'; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    %% Evaluate user choices
    switch act
        case 1
            METHOD = methodarr{nk_input('Select standardization method',0,'m',strjoin(methodarr,'|'),1:numel(methodarr),find(strcmp(methodarr,METHOD)))};
            if ~strcmp(METHOD,'standardization') 
                WINSOPT = []; PX.opt = []; 
            end
        case 2    
            % Define from which data set you compute the mean and std values
            srcfl = nk_input('Compute standardization model using a subgroup of cases?',0,'yes|no',[1,0],sINDDEF);
            if srcfl
                nk_SelectCovariateIndex(NM, NM.TrainParam.FUSION.M, 0);
                sIND = nk_input('Define column indices in the covariate matrix that identify the subgroup(s) for model computation',0,'e',[], Inf );
            else
                sIND = []; CALIBUSE = 2;
            end
        case 3
            % Define to which data set you apply the mean and std values
            dstfl = nk_input('Apply standardization model to a subgroup of cases?',0,'yes|no',[1,0],dINDDEF);
            if dstfl
                if isempty(sIND), n=1; else, n = size(sIND,1); end
                nk_SelectCovariateIndex(NM, NM.TrainParam.FUSION.M, 0);
                if n == 1
                    dIND = nk_input('Define column index in the covariates indicating subgroups to be standardized',0,'e',[], 1);
                else
                    dIND = nk_input('Define column indices in the covariates indicating subgroups to be standardized',0,'e',[], n );
                end
            else
                dIND = [];
            end
        case 4
            % Define whether standardized data should be winsorized
            winsflag = nk_input('Winsorize data (clamping of outliers)?',0,'yes|no',[1,0],WINSOPTDEF);
            if winsflag
                WINSOPT = nk_input('Winsorization +/- threshold (z value)',0,'e',WINSOPT);
                PX = nk_AddParam(WINSOPT, 'winsopt', 1, []);
            else
                WINSOPT = [];PX = nk_AddParam(WINSOPT, 'winsopt', 1, [],'reset');
            end
        case 5
            if IQRFUN == 1, IQRFUN = 2; elseif IQRFUN == 2, IQRFUN = 1; end
        case 6
             if ZeroOut == 1, ZeroOut = 2; elseif ZeroOut == 2, ZeroOut = 1; end
        case 1000
            CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE); 
    end
else
    act = 0;
end

%% Save params
STANDARD.METHOD = METHOD;
STANDARD.sIND = sIND;
STANDARD.dIND = dIND;
STANDARD.WINSOPT = WINSOPT;  
STANDARD.IQRFUN = IQRFUN;
STANDARD.zerooutflag = ZeroOut;
STANDARD.CALIBUSE = CALIBUSE;

% Generate parameter array for preprocessing pipeline runner
if exist('PX','var') && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end
STANDARD.PX = PX;
