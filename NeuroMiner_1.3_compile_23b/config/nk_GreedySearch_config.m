function [ act, param ] = nk_GreedySearch_config(param, SVM, MULTI, defaultsfl, parentstr)

% Minimum # of features
GreedySearch.Direction           = 1;
GreedySearch.EarlyStop.Thresh    = 50;   % Early stopping criterion
GreedySearch.EarlyStop.Perc      = 1;    % Early stopping mode = percentage
GreedySearch.WeightSort          = 1;
GreedySearch.CritGap.flag        = 2;
GreedySearch.CritGap.crit        = 10;

% Wrapper feature evaluation mode at each optimization step
GreedySearch.FeatStepPerc           = 0;
GreedySearch.FeatRandPerc           = 0;
GreedySearch.KneePointDetection     = 2;
GreedySearch.MultiClassOptimization = 1;
GreedySearch.PreSort                = 1;
GreedySearch.PERM.flag              = 0;
GreedySearch.PERM.nperms            = 100;


if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end
if ~defaultsfl
    
    if isfield(param,'GreedySearch') 
        if isfield(param.GreedySearch,'Direction') && ~isempty( param.GreedySearch.Direction)
            GreedySearch.Direction              = param.GreedySearch.Direction;
        end
        if isfield(param.GreedySearch,'WeightSort') && ~isempty( param.GreedySearch.WeightSort)
            GreedySearch.WeightSort             = param.GreedySearch.WeightSort;
        end
        if isfield(param.GreedySearch.EarlyStop,'Thresh') && ~isempty( param.GreedySearch.EarlyStop )
            GreedySearch.EarlyStop.Thresh       = param.GreedySearch.EarlyStop.Thresh; 
            GreedySearch.EarlyStop.Perc         = param.GreedySearch.EarlyStop.Perc;
        end
        if isfield(param.GreedySearch,'FeatStepPerc') && ~isempty( param.GreedySearch.FeatStepPerc)
            GreedySearch.FeatStepPerc           = param.GreedySearch.FeatStepPerc;
        end
        if isfield(param.GreedySearch,'FeatRandPerc') && ~isempty( param.GreedySearch.FeatRandPerc)
            GreedySearch.FeatRandPerc           = param.GreedySearch.FeatRandPerc;
        end
        if isfield(param.GreedySearch,'KneePointDetection') && ~isempty( param.GreedySearch.KneePointDetection)
            GreedySearch.KneePointDetection     = param.GreedySearch.KneePointDetection;
        end
        if isfield(param.GreedySearch,'MultiClassOptimization') && ~isempty( param.GreedySearch.MultiClassOptimization)
            GreedySearch.MultiClassOptimization = param.GreedySearch.MultiClassOptimization;
        end
        if isfield(param.GreedySearch,'PreSort') && ~isempty( param.GreedySearch.PreSort)
            GreedySearch.PreSort = param.GreedySearch.PreSort;
        end
        if isfield(param.GreedySearch,'PERM')
            GreedySearch.PERM                   = param.GreedySearch.PERM;
        end
        if isfield(param.GreedySearch,'CritGap') && ~isempty( param.GreedySearch.CritGap)
            GreedySearch.CritGap = param.GreedySearch.CritGap;
        end
    end
    
    if GreedySearch.Direction == 1
        dirstr = 'Forward';
    else
        dirstr = 'Backward';
    end
    
    % Linear machine?
    if any(strcmp(SVM.kernel.kernstr,{' -t 0', 'lin', 'linear', 'Linear'})) 
        LinMode = 1;
    else
        LinMode = 0;
    end
    
    if LinMode && GreedySearch.MultiClassOptimization ~= 1
        switch GreedySearch.WeightSort
            case 1 
                linstr2 = 'Optimize features for prediction performance';
            case 2
                linstr2 = 'Optimize features for their weights (Guyon''s)';
        end
        linstr = sprintf('|Feature sorting criterion (Linear machine) [ %s ]', linstr2 );
        actind = [1 2];
    else
        linstr = ''; 
        actind = 1;
    end
    
    if GreedySearch.EarlyStop.Thresh
        if GreedySearch.EarlyStop.Perc == 1
            if GreedySearch.Direction == 1
                earlystr = sprintf('Stop when %g%% of the candidate features are still in the pool',GreedySearch.EarlyStop.Thresh);
            else    
                earlystr = sprintf('Stop when %g%% of the features are still used by the model',GreedySearch.EarlyStop.Thresh);
            end
        else
            earlystr = sprintf('Stop at %g selected features',GreedySearch.EarlyStop.Thresh);
        end
    else
        earlystr = 'Early stopping disabled';
    end
    actind = [actind 3];
    randstr = '';
    if GreedySearch.FeatStepPerc
        if GreedySearch.EarlyStop.Perc == 1
            percstr = '%';
        else
            percstr = '';
        end
        stepstr = sprintf('%g%s of top features at each cycle',GreedySearch.FeatStepPerc, percstr);
        if GreedySearch.WeightSort == 1
            if ~GreedySearch.FeatRandPerc
                randoptstr = 'deactivated';
            else
                randoptstr = sprintf('randomly select %g%s  of top features', GreedySearch.FeatRandPerc, percstr);
            end
            randstr = sprintf('|Randomly select features from top-ranked feature block [ %s ]', randoptstr);
            actind = [actind 4 5];
        else
            actind = [actind 4];
        end
    else
        stepstr = 'Each feature will be evaluated';
        actind = [actind 4];
    end
    
    if (~GreedySearch.FeatStepPerc || ~GreedySearch.EarlyStop.Thresh) || ...
         ( GreedySearch.FeatStepPerc > 0 && GreedySearch.EarlyStop.Perc == 1 && (100-GreedySearch.EarlyStop.Thresh) / GreedySearch.FeatStepPerc >= 3 )
        kneestrdef = {'enabled','disabled'};
        kneestr = sprintf('|Kneepoint-based threshold detection [ %s ]',kneestrdef{GreedySearch.KneePointDetection});
        actind = [ actind 6 ];
    else
        kneestr = [];
    end
    
    multistr = []; multipermflagstr = []; multinpermsstr = []; presortstr=[];
    if MULTI.flag && MULTI.train
        if GreedySearch.MultiClassOptimization
            multistrdef = 'Multi-class performance';
        else
            multistrdef = 'Binary/Regression performance';
        end
        multistr = sprintf('|Optimization criterion [ %s ]', multistrdef); actind = [ actind 7 ]; 
        if GreedySearch.MultiClassOptimization
            if GreedySearch.PERM.flag
                multipermflagstrdef = 'activated';
            else
                multipermflagstrdef = 'not activated';
            end
            multipermflagstr = sprintf('|Permutation-based multi-class optimization [ %s ]',multipermflagstrdef ); actind = [actind 8 ];
            if GreedySearch.PERM.flag
                multinpermsstr = sprintf('|Number of permutations to be performed [ %g ]', GreedySearch.PERM.nperms); actind = [ actind 9 ];
            end
            if GreedySearch.PreSort
                presortstrdef = 'yes';
            else
                presortstrdef = 'no';
            end
            presortstr = sprintf('|Sort features according to binary performance before multi-class optimization [ %s ]',presortstrdef); actind = [ actind 10 ];
        end
    end

    critgapstr = [];
    critvaluestr=[];
    if param.datamode == 4
        critgapstrdef = {'enabled','disabled'};
        critgapstr = sprintf('|Critical-gap based early stopping [ %s ]',critgapstrdef{GreedySearch.CritGap.flag}); actind = [ actind 11 ];
        if GreedySearch.CritGap.flag == 1
            critvaluestr = sprintf('|Early stop-inducing gap between CV1 training and test criteria in %% [ %g ]', GreedySearch.CritGap.crit); actind = [ actind 12 ];
        end
    end
    
    mestr = 'Greedy feature selection setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
    
    nk_PrintLogo
    act = nk_input('Define Greedy Search parameters',0,'mq', ...
                    [sprintf('Search direction [ %s ]', dirstr) ...
                     linstr ...
                     sprintf('|Early stopping [ %s ]', earlystr) ...
                     sprintf('|Feature stepping [ %s ]', stepstr) ...
                     randstr ...
                     kneestr ...
                     multistr ...
                     multipermflagstr ...
                     multinpermsstr ...
                     presortstr ...
                     critgapstr ...
                     critvaluestr], actind);
    switch act
        case 1
            GreedySearch.Direction = nk_input('Feature search mode',0,'m', ...
                'Forward selection|Backward selection',[1,2], GreedySearch.Direction);
        case 2
            GreedySearch.WeightSort = nk_input('Feature sorting criterion',0,'m', ...
                'Performance|Weights (Guyon''s method)',[1,2], GreedySearch.WeightSort);
            if GreedySearch.WeightSort == 2, GreedySearch.FeatRandPerc = 0; end
        case 3
            GreedySearch.EarlyStop.Thresh = nk_input('Early stopping threshold (0 = disables early stopping)', 0,'e', ...
                GreedySearch.EarlyStop.Thresh);
            if GreedySearch.EarlyStop.Thresh
                GreedySearch.EarlyStop.Perc = nk_input('Absolute number or Percentage of features',0, 'Percentage|Absolute', ...
                    [1,2], GreedySearch.EarlyStop.Perc);
            end
        case 4
            GreedySearch.FeatStepPerc = nk_input('Define % of top-ranked features selected at each wrapper cycle (0 = disables block selection)',0,'e',GreedySearch.FeatStepPerc);
        case 5
            GreedySearch.FeatRandPerc = nk_input('Randomly select % of features in block of top-ranked features (0 = disables random selection)',0,'e',GreedySearch.FeatRandPerc);
        case 6
            if GreedySearch.KneePointDetection == 1, GreedySearch.KneePointDetection = 2; else, GreedySearch.KneePointDetection = 1; end
        case 7
            GreedySearch.MultiClassOptimization = ~GreedySearch.MultiClassOptimization;
        case 8
            GreedySearch.PERM.flag = ~GreedySearch.PERM.flag;
        case 9
            GreedySearch.PERM.nperms = nk_input('Define number of permutations for multi-class feature optimization',0,'i', GreedySearch.PERM.nperms);     
        case 10
            GreedySearch.PreSort = ~GreedySearch.PreSort; 
        case 11
            if GreedySearch.CritGap.flag == 1, GreedySearch.CritGap.flag = 2; else, GreedySearch.CritGap.flag = 1; end
        case 12
            GreedySearch.CritGap.crit = nk_input('Define critical gap in %% between CV1 training and test data for early stopping',0,'i', GreedySearch.CritGap.crit);     
    end
else
    act = 0;
end

param.GreedySearch = GreedySearch; 
