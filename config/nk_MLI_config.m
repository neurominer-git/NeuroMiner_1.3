% =========================================================================
% FORMAT param = nk_MLI_config(param)
% =========================================================================
% NeuroMiner menu configurator for the interpretation of model predictions 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2022
function [MLI, act] = nk_MLI_config(MLI, M, defaultsfl, parentstr)
global NM

if ~exist("defaultsfl","var") || isempty(defaultsfl)
    defaultsfl=false;
end

method          = 'posneg';
upper_thresh    = 95;
lower_thresh    = 5;
nperms          = 1024; % 2^10, i.e. for up to 10 features all coalitions will be covered in shapley
max_iter        = 1000;
n_visited       = 100;
frac            = .1;
samplefrac      = .25;
usemap          = 0;
mapfeat         = 'cvr';
cutoff          = [-2 2];
cutoffmode      = 'absolute';
cutoffoperator  = 1;
znormdata       = 1;
isimaging       = false;
MLIatlasflag    = 0;
IO = NM.datadescriptor{M}.input_settings;
if ~strcmp(NM.datadescriptor{M}.source,'matrix')
    isimaging   = true;
end
if isimaging
    MLIatlasflag = 1;
    MLIatlasfile = [];
    MLIatlasvec  = [];
    MLIcsvfile   = [];
    MLIcsvdelim  = 'space';
    MLIcsvcol    = 'ROIdesc';
    MLIcsvid     = 'ROIid';
    MLIcsvdesc   = [];
    MLIcsvnum    = [];
end
% samples = 'All';
% sampleID = '';

if ~defaultsfl

    if ~exist('MLI','var') || isempty(MLI)
        MLI = nk_MLI_config([], M, 1); 
    end
    
    method          = MLI.method;
    nperms          = MLI.nperms;
    upper_thresh    = MLI.upper_thresh;
    lower_thresh    = MLI.lower_thresh;
    max_iter        = MLI.max_iter;
    n_visited       = MLI.n_visited;
    znormdata       = MLI.znormdata;
    if numel(MLI.Modality) >= M && ~isempty(MLI.Modality{M}) 
        frac            = MLI.Modality{M}.frac;
        usemap          = MLI.Modality{M}.MAP.flag;
        mapfeat         = MLI.Modality{M}.MAP.map;
        cutoff          = MLI.Modality{M}.MAP.cutoff;
        cutoffmode      = MLI.Modality{M}.MAP.percentmode;
        cutoffoperator  = MLI.Modality{M}.MAP.operator;
        if isimaging
            MLIatlasflag = MLI.Modality{M}.imgops.flag;
            MLIatlasfile = MLI.Modality{M}.imgops.atlasfile;
            MLIatlasvec  = MLI.Modality{M}.imgops.atlasvec;
            MLIcsvfile   = MLI.Modality{M}.imgops.csvfile;
            MLIcsvdelim  = MLI.Modality{M}.imgops.csvdelim;
            MLIcsvid     = MLI.Modality{M}.imgops.csvid ;
            MLIcsvcol    = MLI.Modality{M}.imgops.csvcol;
            MLIcsvdesc   = MLI.Modality{M}.imgops.csvdesc;
            MLIcsvnum    = MLI.Modality{M}.imgops.csvnum;        
        end

        % if isfield(MLI,'samples')
        %     samples = MLI.samples;
        %     sampleID = MLI.sampleID;
        % end
        if isfield(MLI,'frac')
            samplefrac = MLI.frac;
        end
    end

    switch method
        case 'posneg'
            MethodStr = sprintf('Upper/lower percentiles (%g%%/%g%%)', upper_thresh, lower_thresh);
            OcclusionUpperThreshStr = sprintf('|Define upper percentile [ %g ]', upper_thresh);
            OcclusionLowerThreshStr = sprintf('|Define lower percentile [ %g ]', lower_thresh);
            DEFMETHOD = 1;
            % DEFSAMPLES = 1;
            mnuact = [1 2 3];
        case 'median'
            MethodStr = 'Median';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 2;
            % DEFSAMPLES = 1;
            mnuact = 1;
        case 'medianflip'
            MethodStr = 'Median span flipped';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 3;
            % DEFSAMPLES = 1;
            mnuact = 1;
        case 'medianmirror'
            MethodStr = '100-quantile (median mirror)';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 4;
            %DEFSAMPLES = 1;
            mnuact = 1;
        case 'random'
            MethodStr = 'Random value';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 5;
            % DEFSAMPLES = 1;
            mnuact = 1;
        case 'shapley'
            MethodStr = 'Shapley values';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 6;
            % DEFSAMPLES = 1;
            mnuact = 1;
%         case 'tree'
%             MethodStr = 'Treeinterpreter';
%             OcclusionUpperThreshStr = '';
%             OcclusionLowerThreshStr = '';
%             DEFMETHOD = 6;
%             DEFSAMPLES = 1;
%             mnuact = 1;
    end
    OcclusionMethodStr = ['Define occlusion method [ ' MethodStr ' ]']; 
    
    if isinf(nperms)
        IterStr = 'Automated stopping';
    else
        IterStr = sprintf('%g iterations', nperms);
    end
    OcclusionIterStr = ['|Define no. of iterations [ ' IterStr ' ]']; 
    mnuact = [mnuact 4];
    
    if isinf(nperms) && ~isequal(method, 'shapley')
        OcclusionVisitedStr = ['|Define minimum number of feature visits [ ' num2str(n_visited) ' ]'];
        mnuact = [ mnuact 5 ];
    else
        OcclusionVisitedStr = '';
    end
    if ~isequal(method, 'shapley')
        OcclusionFracStr = ['|Define fraction (0-1) of features to be visited per iteration [ ' num2str(frac) ' ]'];
        mnuact = [ mnuact 6 ];
    else
        OcclusionFracStr = '';
    end
    
    if ~usemap 
        MapFlagStr = 'no map specified'; 
    else
        MapFlagStr = 'yes';
    end
    OcclusionMapFlagStr = ['|Use map from model visualization to operate in pre-determined feature space [ ' MapFlagStr ' ]'];
    mnuact = [ mnuact 7 ];
    if usemap
        switch mapfeat
            case 'cvr'
               MapFeatStr = 'Cross-validation ratio map'; DEFMAP = 1; 
            case 'p_sgn'
               MapFeatStr = 'Sign-based consistency map (p value, uncorrected)'; DEFMAP = 2;
            case 'p_FDR_sgn'
               MapFeatStr = 'Sign-based consistency map (p value, FDR-corrected)'; DEFMAP = 3;   
        end
        OcclusionMapFeatStr = ['|Define which map to use [ ' MapFeatStr ' ]'];
        mnuact = [ mnuact 8 ];
        OcclusionMapCutoffStr = ['|Define map cutoff value (scalar or 2-value vector) [ ' num2str(cutoff) ' ]'];
        mnuact = [ mnuact 9 ];
        switch cutoffmode
            case 'absolute'
                DEFCUTOFFMODE = 1;
            case 'percentile'
                DEFCUTOFFMODE = 2;
        end
        OcclusionMapCutoffModeStr = ['|Define map cutoff method [ ' cutoffmode ' ]'];
        mnuact = [ mnuact 10 ];
        if numel(cutoff)>1
            MapCutoffOperatorStr = {'<,>', '<=,>=', '>,<', '>=,<='};
        else
            MapCutoffOperatorStr = {'<', '<=', '>', '>=', '==', '~='};
        end
        OcclusionMapCutoffOperator = ['|Define map cutoff operator [ ' MapCutoffOperatorStr{cutoffoperator} ' ]'];
        mnuact = [ mnuact 11 ];
    else
        OcclusionMapFeatStr = '';
        OcclusionMapCutoffStr = '';
        OcclusionMapCutoffModeStr = '';
        OcclusionMapCutoffOperator = '';
    end

    if ~isequal(method, 'shapley')
        ZnormDataOpts = {'No, only group-level normalization', ...
            'Yes, mean centering at the case level', ...
            'Yes, Z-normalization at the case level', ...
            'Yes, scaling to [-1,1] at the case level'};
        OcclusionZnormData = ['|Normalize prediction change estimates at the case level [ ' ZnormDataOpts{znormdata} ' ]' ];
        mnuact = [ mnuact 12 ];
    else
        OcclusionZnormData = '';
    end
    OcclusionAtlasFlag = '';
    OcclusionAtlasImgFile = '';
    OcclusionAtlasCSVDelim = '';
    OcclusionAtlasCSVID = '';
    OcclusionAtlasCSVCol = '';
    OcclusionAtlasCSVFile = '';
    OcclusionAtlasCSVFileRead = '';
    if isimaging
        if ~MLIatlasflag
            OcclusionAtlasFlag = '|Bind data modification to neuroanatomical atlas [ no ]' ; 
            mnuact = [ mnuact 13 ];
        else
            if exist(MLIatlasfile,"file")
                [~,AtlasFilename, Filename_ext] = fileparts(MLIatlasfile);
            else
                AtlasFilename = 'undefined'; Filename_ext  = '';
            end
            if exist(MLIcsvfile,"file")
                [~,AtlasCSVname, CSVname_ext] = fileparts(MLIcsvfile);
            else
                AtlasCSVname = 'undefined'; CSVname_ext = '';
            end
            OcclusionAtlasFlag =    '|Bind data modification to neuroanatomical atlas [ yes ]' ;
            OcclusionAtlasImgFile =  [ '|Define & import atlas imaging file [ ' AtlasFilename Filename_ext ' ]' ];
            OcclusionAtlasCSVFile =  [ '|Define atlas descriptor file [ ' AtlasCSVname CSVname_ext ' ]' ];
            OcclusionAtlasCSVDelim = [ '|Define delimiter used in atlas descriptor file [ ' char(MLIcsvdelim) ' ]' ];
            OcclusionAtlasCSVID =    [ '|Define ID column in atlas descriptor file [ ' MLIcsvid ' ]' ];
            OcclusionAtlasCSVCol =   [ '|Define descriptor column in atlas descriptor file [ ' MLIcsvcol ' ]' ];
            if ~strcmp(MLIcsvfile,'undefined')
                OcclusionAtlasCSVFileRead = '|Import atlas descriptor file';
                mnuact = [ mnuact 13:19 ];
            else
                mnuact = [ mnuact 13:18 ];
            end
            
        end
    end

    % OcclusionSamplesStr = ['|Define for which sample(s) to interpret the model [ ' samples ' ]|'];
    % mnuact = [ mnuact 20 ];
    % 
    % if isequal(samples, 'User-defined')
    %     if isempty(sampleID)
    %         sampleIDstr = '';
    %     elseif numel(sampleID) > 4
    %         sampleIDstr = [num2str(numel(sampleID)) ' IDs'];
    %     else
    %         sampleIDstr = strjoin(sampleID);
    %     end
    %     OcclusionSampleIDStr = ['|Specify sample ID(s) for which to interpret the model [ ' sampleIDstr ' ]|'];
    %     mnuact = [ mnuact 21 ];
    % else
    %     OcclusionSampleIDStr = '';
    % end


    if isequal(method,'shapley')
        SampleFracStr = ['|Define fraction (0-1) of training samples to be visited per iteration [ ' num2str(samplefrac) ' ]'];
        mnuact = [mnuact 22];
    else
        SampleFracStr = '';
    end
    
    mnustr = [ OcclusionMethodStr ...
              OcclusionUpperThreshStr ...
              OcclusionLowerThreshStr ...
              OcclusionIterStr ...
              OcclusionVisitedStr ...
              OcclusionFracStr ...
              OcclusionMapFlagStr ...
              OcclusionMapFeatStr ...
              OcclusionMapCutoffStr ...
              OcclusionMapCutoffModeStr ...
              OcclusionMapCutoffOperator ...
              OcclusionZnormData ...
              OcclusionAtlasFlag ...
              OcclusionAtlasImgFile ...
              OcclusionAtlasCSVFile ...
              OcclusionAtlasCSVDelim ...
              OcclusionAtlasCSVID ...
              OcclusionAtlasCSVCol ...
              OcclusionAtlasCSVFileRead ...
              SampleFracStr ];
              % OcclusionSamplesStr ...
              % OcclusionSampleIDStr ...
              
    
    nk_PrintLogo
    if numel(MLI.Modality) >= M && isfield(MLI.Modality{M},'imgops') && ~isempty(MLI.Modality{M}.imgops) && ...
            MLI.Modality{M}.imgops.flag && ~isempty(MLIcsvfile) && isfield(MLI.Modality{M}.imgops,'T') 
        if isfield(MLI.Modality{M}.imgops,'csvnum') && ~isempty(MLI.Modality{M}.imgops.csvnum)
            fprintf('\nFirst line of atlas descriptor file [imported]: %s\n\n', MLIcsvfile); 
        else
            fprintf('\nFirst line of atlas descriptor file: %s\n\n', MLIcsvfile); 
        end
        disp(MLI.Modality{M}.imgops.T(1,:));
    end
    if M>1
        mestr = sprintf('Model interpretation parameters [ Modality #%g: %s ]', M, NM.datadescriptor{M}.desc) ; 
    else
        mestr = 'Model interpretation parameters'; 
    end
    navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', mnustr, mnuact);
    
    switch act
        case 1
            method = char(nk_input('Define occlusion method', 0, 'm', ...
                ['Replace feature with upper and lower percentile cutoffs|' ...
                 'Replace feature with median|' ...
                 'Replace feature by adding/subtracting the median quantile|' ...
                 'Replace feature by its mirror value at the median quantile (100-quantile approach)|' ...
                 'Replace feature by randomly picking a value from the feature''s training sample distribution|' ...
                 'Replace features according to Shapley'], ...
                 {'posneg','median','medianflip','medianmirror','random','shapley'}, DEFMETHOD));
%                  'Replace features according to Shapley|' ...
%                  'Use treeinterpreter'], ...
%                  {'posneg','median','medianflip','random','shapley','tree'}, DEFMETHOD));
        case 2
            upper_thresh = nk_input('Define upper percentile', 0, 'e', upper_thresh);
        case 3
            lower_thresh = nk_input('Define upper percentile', 0, 'e', lower_thresh);
        case 4
            nperms = nk_input('Define no. of iterations [inf for automated mode]', 0, 'e', nperms);
        case 5
            n_visited = nk_input('Define minimum number of visits per feature for automated mode', 0, 'i', n_visited);
        case 6
            frac = nk_input('Define fraction of data space to be visited per iteration', 0, 'e', frac);
        case 7
            if ~usemap, usemap=1; else, usemap = 0; end
        case 8
            mapfeat = char(nk_input('Define occlusion method', 0, 'm', ...
                ['Cross-validation ratio map|' ...
                 'Sign-based consistency map (p value, uncorrected)|' ...
                 'Sign-based consistency map (p value, FDR-corrected)'], ...
                 {'cvr','p_sgn','p_FDR_sgn'}, DEFMAP));
        case 9
            switch mapfeat
                case 'cvr'
                   MapFeatStr = 'Cross-validation ratio map';
                case 'p_sgn'
                   MapFeatStr = 'Sign-based consistency map (p value, uncorrected)';
                case 'p_FDR_sgn'
                   MapFeatStr = 'Sign-based consistency map (p value, FDR-corrected)';
            end
            cutoff = nk_input(['Define cutoff(s) to be applied to ' MapFeatStr ], 0, 'e', cutoff);
        case 10
            cutoffmode = char(nk_input('Define map thresholding mode', 0, 'm', ...
                ['Absolute mode|' ...
                 'Percentile mode' ], ...
                 {'absolute','percentile'}, DEFCUTOFFMODE));
        case 11
            if numel(cutoff) > 1

                cutoffoperator = nk_input('Define map thresholding operator', 0, 'm', ...
                    ['<,>|' ...
                     '<=,>=|' ...
                     '>,<|' ...
                     '>=,<='], ...
                     1:4, cutoffoperator);
            else
                cutoffoperator = nk_input('Define map thresholding operator', 0, 'm', ...
                    ['<|' ...
                     '<=|' ...
                     '>|' ...
                     '>=|' ...
                     '==|', ...
                     '~='], ... 
                     1:6, cutoffoperator);
            end
        case 12
            znormdata = nk_input('Define centering procedure', 0, 'm', 'None|Mean centering|Z-normalization|Scaling to [-1, 1]', 1:4, znormdata);
        case 13
            MLIatlasflag = ~MLIatlasflag;
        case 14
            [MLIatlasfile, Vol] = nk_FileSelector(1, IO.datasource, 'Select atlas file', IO.filt, MLIatlasfile);
            imgtype = NM.datadescriptor{M}.input_settings.datasource;
            switch imgtype
                case {'spm','nifti'}
                     Thresh = NM.datadescriptor{M}.input_settings.Thresh;
                     Vm = spm_vol(NM.brainmask{M}); 
                     MLIatlasvec = round(nk_ReturnSubSpaces(Vol, Vm, 1, 1, Thresh),0);
                case 'surf'
                     MLIatlasvec = round(MRIread(MLIatlasfile),0);
            end
        case 15
            MLIcsvfile = nk_FileSelector(1,0,'Select atlas descriptor file', '.*\.txt$|.*\.dat$|.*\.csv', MLIcsvfile);
            try
                 MLI.Modality{M}.imgops.T = readtable(MLIcsvfile, 'Delimiter', MLIcsvdelim);
            catch ERR
                errordlg(sprintf('The atlas descriptor file %s could not be read!\nError: %s\nCheck your settings!', ...
                    MLIcsvfile, ERR.message), 'MLI configuration module');
            end
        case 16
            MLIcsvdelimopts = {';',',','tab','space'};
            MLIcsvdelim = spm_input('Choose the delimiter', 0, 'm', strjoin(MLIcsvdelimopts,'|'), ...
                MLIcsvdelimopts, find(strcmp(MLIcsvdelimopts,MLIcsvdelim)));
            try
                 MLI.Modality{M}.imgops.T = readtable(MLIcsvfile, 'Delimiter', MLIcsvdelim);
            catch ERR
                errordlg(sprintf('The atlas descriptor file %s could not be read!\nError: %s\nCheck your settings!', ...
                    MLIcsvfile, ERR.message), 'MLI configuration module');
            end
        case 17
            MLIcsvid = spm_input('Define the name of the table column containing the ROI indices', 0,'s', MLIcsvid);
        case 18
            MLIcsvcol = spm_input('Define the name of the table column containing the ROI labels', 0,'s', MLIcsvcol);
        case 19
            try
                MLIcsvdesc = MLI.Modality{M}.imgops.T.(MLIcsvcol);
                MLIcsvnum =  MLI.Modality{M}.imgops.T.(MLIcsvid);
            catch ERR
                errordlg(sprintf('The atlas descriptor file %s could not be processed!\nError: %s\nCheck your settings!', ...
                    MLIcsvfile, ERR.message), 'MLI configuration module');
            end
        % case 20
        %     samples = char(nk_input('Define samples as an array of IDs', 0, 'm', ...
        %         ['All|' ...
        %          'User-defined' ], ...
        %          {'All','User-defined'}, DEFSAMPLES));
        % case 21
        %     try 
        %         sampleID = strjoin(sampleID);
        %     catch
        %         sampleID = sampleID;
        %     end
        %     while ~iscellstr(sampleID)
        %         sampleID = nk_input('Enter sample IDs as a cell array of strings (e.g. {''ID1'',''ID2''})', 0, 'e', sampleID);
        %     end
        %     %make sure the cell array of strings is in one column
        %     sampleID = reshape(sampleID,size(sampleID,1)*size(sampleID,2),1);

        case 22
            samplefrac = nk_input('Define fraction of training samples to be visited per iteration', 0, 'e', samplefrac);
    end
end
MLI.method = method;
MLI.nperms = nperms;
MLI.upper_thresh = upper_thresh;
MLI.lower_thresh = lower_thresh;
MLI.max_iter = max_iter;
MLI.n_visited = n_visited;
MLI.znormdata = znormdata;
MLI.frac = samplefrac;
MLI.Modality{M}.frac = frac;
MLI.Modality{M}.MAP.flag = usemap;
MLI.Modality{M}.MAP.map = mapfeat;
MLI.Modality{M}.MAP.cutoff = cutoff;
MLI.Modality{M}.MAP.percentmode = cutoffmode;
MLI.Modality{M}.MAP.operator = cutoffoperator;
if isimaging
    MLI.Modality{M}.imgops.flag = MLIatlasflag;
    MLI.Modality{M}.imgops.atlasfile = MLIatlasfile;
    MLI.Modality{M}.imgops.atlasvec = MLIatlasvec;
    MLI.Modality{M}.imgops.atlasvec_unique = unique(MLIatlasvec);
    MLI.Modality{M}.imgops.csvfile = MLIcsvfile;
    MLI.Modality{M}.imgops.csvdelim = MLIcsvdelim;
    MLI.Modality{M}.imgops.csvid = MLIcsvid;
    MLI.Modality{M}.imgops.csvcol = MLIcsvcol;
    MLI.Modality{M}.imgops.csvdesc = MLIcsvdesc;
    MLI.Modality{M}.imgops.csvnum = MLIcsvnum;
else
    MLI.Modality{M}.imgops = [];
end
% MLI.samples = samples;
% try
%     MLI.sampleID = sampleID_cell;
% catch
%     MLI.sampleID = sampleID;
% end
