function [DR, PX, TEMPLPROC, USEALL, act] = nk_DimRed_main_config(DR, PX, TEMPLPROC, USEALL, parentstr, defaultsfl)
global EXPERT 

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
CALIBUSE = 2; 

if ~defaultsfl
    if isempty(DR), [DR, PX, TEMPLPROC, USEALL] = nk_DimRed_main_config(DR, PX, TEMPLPROC, USEALL, parentstr, true); end
    if isfield(DR,'RedMode'),   RedMode = DR.RedMode; else, RedMode = []; end
    if isfield(DR,'dims'),      dims = DR.dims; else, dims = []; end
    if isfield(DR,'PercMode'),  PercMode = DR.PercMode; else, PercMode = []; end
    if isfield(DR,'CALIBUSE'),  CALIBUSE = DR.CALIBUSE; end
    if ~exist('TEMPLPROC','var') || isempty(TEMPLPROC), TEMPLPROC = false; end
    if ~exist('USEALL','var') || isempty(USEALL), USEALL = false; end
    if ~ischar(RedMode)
        DRSTR_REDMODE = 'undefined';
    else
        DRSTR_REDMODE = RedMode;
    end
    
    if strcmp(DRSTR_REDMODE,'PCA')
        if isempty(PercMode)
            DRSTR_PERCMODE = '(undefined)';
        else
            switch PercMode
                case 1
                    DRSTR_PERCMODE = ' (Absolute)';

                case 2
                    DRSTR_PERCMODE = ' (Percentage)';

                case 3
                    DRSTR_PERCMODE = ' (Energy)';
            end
        end
    else
        DRSTR_PERCMODE = '';
    end
    
    if isempty(dims)
        DRSTR_DIM = 'undefined';
    else
        DRSTR_DIM = [ nk_ConcatParamstr(dims) DRSTR_PERCMODE ];
    end

    if strcmp(DRSTR_REDMODE,'PLS') 
        menustr = ['Select dimensionality reduction method [ ' DRSTR_REDMODE ' ]'];
        menuact = 1;
    elseif strcmp(DRSTR_REDMODE, 'fastICA')
        menustr = ['Select dimensionality reduction method [ ' DRSTR_REDMODE ' ]|', ...
                   'Define fastICA settings [ dimensions: ' DRSTR_DIM ' ]'];
        menuact = [1 3];
    else 
        menustr = ['Select dimensionality reduction method [ ' DRSTR_REDMODE ' ]|', ...
                   'Define dimenisonality of mapping results [ ' DRSTR_DIM ' ]'];
        menuact = 1:2;
    end
    % Configuration of matrix realignment only available for expert
    % settings
    if EXPERT && contains(DRSTR_REDMODE,{'PCA', 'RobPCA', 'FactorAnalysis', 'SparsePCA', 'ProbPCA'})
        if TEMPLPROC, DRSTR_PROCRUST = 'yes'; else, DRSTR_PROCRUST = 'no'; end 
        menustr = [menustr '|Use Procrustes rotation to align CV1-level factorizations [ ' DRSTR_PROCRUST ' ]'];
        menuact = [menuact 4];
        if USEALL==1, DRSTR_PROCRUST_USEALL = 'yes'; else, DRSTR_PROCRUST_USEALL = 'no'; end
        menustr = [menustr '|Use all the training data for template generation [ ' DRSTR_PROCRUST_USEALL ' ]'];
        menuact = [menuact 5];
    end
    
    [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);
    
    nk_PrintLogo
    mestr = 'Dimensionality reduction'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    switch act
        case 1
            [DR, PX] = return_dimred(DR, PX, parentstr);
            [DR.dims, DR.PercMode] = nk_ExtDim_config(DR.RedMode, [], [], true);
        case 2 % for compatibility
            t_act = 1; while t_act > 0, [dims, PercMode, t_act] = nk_ExtDim_config(RedMode, PercMode, dims, 0, navistr);end
            DR.dims = dims; DR.PercMode = PercMode;
        case 3
            act = 1;
            while act > 0
                [DR, PX, act] = cv_fastICA_config(DR, PX, parentstr, defaultsfl);
            end
        case 4
            TEMPLPROC = ~TEMPLPROC; 
        case 5
            USEALL = ~USEALL;
        case 1000
            CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE);
    end
else
    [DR, PX] = return_dimred(DR, PX, parentstr, true);
    [DR.dims, DR.PercMode] = nk_ExtDim_config(DR.RedMode, [], [], true);
    act = 0;
end

% Add parameters for dimensionality reduction optimization
PX = nk_AddParam(DR.dims, 'dimensions', 1, PX);
if numel(DR.dims)<2, DR.EXTOPT = true; else, DR.EXTOPT = false; end
DR.CALIBUSE = CALIBUSE;

end

% _________________________________________________________________________

function [DR, PX] = return_dimred(DR, PX, parentstr, defaultsfl)
global EXPERT NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl

   
    menustr = ['Principal Component Analysis                   (PCA)|' ...
               'Robust Principal Component Analysis            (significantly slower than PCA but more accurate in low N problems)|' ...
               'Non-negative Matrix Factorization (NNMF)       (significantly slower than PCA but more parsimonious)|' ...
               'Orthogonal Projective NNMF                     (significantly slower but more robust compared to NNMF)|' ... 
               'Fast-than-fast NNMF                            (fastest NNMF through random prejections and Nesterov iterations)|' ...
               'Partial Least Squares                          (PLS)|' ...
               'Sparse PCA                                     (... a ''statistical segmentation'' tool)|' ...
               'Probabilistic PCA                              (low-D only, >32 GB RAM for high-D data)|' ...
               'Factor analysis                                ()|' ...
               'Locality Preserving Projections                (low-D only, #cases > #features)|' ...
               'Linear Local Tangent Space Alignment           (low-D only, >32 GB RAM for high-D data)|' ...
               'Large-Margin Nearest Neighbour                 ()|' ...
               'Linear Discriminant Analysis                   (low-D only)|' ...
               'Neighborhood Component Analysis                (low-D only)|' ...
               'Independent Component Analysis                  (fastICA)'];
    menuact = { 'PCA', ...
                'RobPCA',... 
                'NNMF', ...
                'optNMF', ...
                'NeNMF', ...
                'PLS', ...
                'SparsePCA', ...
                'ProbPCA', ...
                'FactorAnalysis',...
                'LPP',...
                'LLTSA',...
                'LMNN',...
                'LDA',...
                'NCA',...
                'fastICA'};    
    if isfield(DR,'RedMode'), def = find(strcmp(menuact,DR.RedMode)); else, def = 1; end
    if EXPERT
        menustr = [menustr '|' ...
               '* General Discriminant Analysis                (low-D only)|' ...
               '* Laplacian Eigenmaps                          ()|' ...
               '* Local Linear Embedding                       ()']; 
           menuact = [ menuact, ...
                'GDA',...
                'Laplacian',...
                'LLE'];
    end
    nk_PrintLogo
    mestr = 'Dimensionality reduction techniques'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    if isempty(def), def=1; end
    if EXPERT
        fprintf('\n*  No reconstruction of data in original space supported!');
    end
    act = char(nk_input(mestr,0,'mq', menustr, menuact, def));
    if ~strcmp(act,'BACK'), DR.RedMode = act; end
    if ~exist('PX','var'), PX = []; end

    switch act

        case 'PCA'
            DR.DRsoft = 1; PX = nk_AddParam([], [], [], PX,'reset');
            
        case 'SparsePCA'
            DR.DRsoft = 1;
            if isfield(DR,'SparsePCA')
                lambda = DR.SparsePCA.lambda;
                stop = DR.SparsePCA.stop;
                maxiter = DR.SparsePCA.maxiter;
            else
                lambda = Inf;
                stop = [ .1 .2 .3 .4 .5 ];
                maxiter = 300;
            end
            DR.SparsePCA.lambda     =   nk_input('Lambda parameter (positive ridge term coefficient; set to Inf for high-D problems)',0,'e', lambda);
            DR.SparsePCA.stop       =   nk_input('Stop parameter (if positive => upper bound on the L1-norm of the BETA coefficients, if 0 => regular PCA)',0,'e', stop);
            DR.SparsePCA.maxiter    =   nk_input('max # of iterations)',0,'w1', maxiter);
            PX = nk_AddParam(DR.SparsePCA.lambda, 'SparsePCA-Lambda', 1, PX, 'replace');
            PX = nk_AddParam(DR.SparsePCA.stop, 'SparsePCA-Stop', 1, PX);
            PX = nk_AddParam(DR.SparsePCA.maxiter, 'SparsePCA-Iter', 1, PX);
            
        case {'RobPCA','FactorAnalysis','LDA','GDA','MDS'}
            PX = nk_AddParam([], [], [], PX,'reset');
            
        case 'PLS'
            PX = nk_AddParam([], [], [], PX,'reset');
            if ~isfield(DR,'PLS'), [~, DR.PLS ] = nk_PLS_config(NM, [], [], true); end
            act = 1; while act >0, [act, DR.PLS, PX ] = nk_PLS_config(NM, DR.PLS, PX, false, navistr); end
  
        case 'NNMF'
            
            PX = nk_AddParam([], [], [], PX,'reset');
            NMFmethod               = {'nmf','sparsenmf','convexnmf','orthnmf'};
            NMFmethoddesc           = {'Standard Non-Negative Matrix Factorization',...
                                      'Sparse Non-Negative Matrix Factorization', ...
                                      'Convex Non-Negative Matrix Factorization', ...
                                      'Orthogonal Non-Negative Matrix T-Factorization'};
            if ~isfield(DR,'NMFmethod'), DR.NMFmethod = NMFmethod{1}; end
            DR.NMFmethod            = NMFmethod{nk_input('Choose NMF method',0,'m',strjoin(NMFmethoddesc,'|'),1:numel(NMFmethoddesc), find(strcmp(NMFmethod,DR.NMFmethod)))};
            DR.options.dis          = 0;
            
        case 'NeNMF'
            NMFmethod               = {'RPI_NeNMF','RSI_NeNMF','VANILLA_NeNMF'};
            NMFmethoddesc           = { 'Randomized Power Iterations NeNMF', ...
                                        'Non-negative Matrix Factorization via Nesterov''s Optimal Gradient Method: Improved via Randomized Subspace Iterations NeNMF', ...
                                        'Non-negative Matrix Factorization via Nesterov''s Optimal Gradient Method' };
            if ~isfield(DR,'NMFmethod'), DR.NMFmethod = NMFmethod{1}; end
            if ~isfield(DR,'tmax'), DR.tmax= 60; end
            defitem = find(strcmp(NMFmethod,DR.NMFmethod));
            if isempty(defitem), defitem = NMFmethoddesc{1}; end
            DR.NMFmethod            = NMFmethod{nk_input('Choose NMF method',0,'mq',strjoin(NMFmethoddesc,'|'),1:numel(NMFmethoddesc), defitem)};
            DR.tmax                 = nk_input('Maximum run-time in seconds',0,'i',DR.tmax);
        case 'ProbPCA'
            if isfield(DR,'ProbPCA'), iter= DR.ProbPCA.iter; else, iter = 200; end
            DR.ProbPCA.iter = nk_input('# of iterations ',0,'e',iter);
            PX = nk_AddParam(DR.ProbPCA.iter, 'ProbPCA-Iter', 1, PX, 'replace');

        case 'LPP'
            if isfield(DR,'LPP') 
                Sigma = DR.LPP.Sigma; 
                K = DR.LPP.k; 
            else
                Sigma = 1; K = 12;
            end
            DR.LPP.Sigma = nk_input('Bandwith of the Gaussian kernel (Sigma)',0,'e',Sigma);
            DR.LPP.k = nk_input('# of nearest neighbors (k)',0,'e',K);
            PX = nk_AddParam(DR.LPP.Sigma, 'LPP-Sigma', 1, PX, 'replace');
            PX = nk_AddParam(DR.LPP.k, 'LPP-K', 1, PX);

        case 'LLTSA'
            if isfield(DR,'LLTSA'), k=DR.LLTSA.k; else, k=12; end
            DR.LLTSA.k = nk_input('# of nearest neighbors (k)',0,'e',k);                 
            PX = nk_AddParam(DR.LLTSA.k, 'LLTSA-K', 1, PX, 'replace');

        case 'Laplacian'
            if isfield(DR,'Laplacian'), k=DR.Laplacian.k; else, k=12; end
            DR.Laplacian.k = nk_input('# of nearest neighbors (k)',0,'e',k);                 
            PX = nk_AddParam(DR.Laplacian.k, 'Laplacian-K', 1, PX, 'replace');
            DR.Laplacian.s = nk_input('Sigma',0,'e',1);              
            PX = nk_AddParam(DR.Laplacian.s, 'Laplacian-Sigma', 1, PX);

        case 'LLE'
            if isfield(DR,'LLE'), k = DR.LLE.k; else, k=12; end
            DR.LLE.k = nk_input('# of nearest neighbors (k)',0,'e',k);
            PX = nk_AddParam(DR.LLE.k, 'LLE-K', 1, PX, 'replace');
        
        case 'NCA'
            DR.NCA.lambda = nk_input('Lambda range',0,'e',0.5);
            PX = nk_AddParam(DR.NCA.lambda, 'NCA-Lambda', 1, PX, 'replace');

        case 'LMNN'
            PX = nk_AddParam([], [], [], PX,'reset');

        case 'fastICA'
            DR.DRsoft = 1; PX = nk_AddParam([], [], [], PX,'reset');
       
    end
else
    DR.RedMode = 'PCA'; DR.DRsoft = 1; PX = nk_AddParam([], [], [], PX,'reset');
end

end