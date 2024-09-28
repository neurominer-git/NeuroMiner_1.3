function [act, param] = nk_PFA_config(param, defaultsfl, parentstr)

PFA.N                = 10;   % Number of particles
PFA.max_Iter         = 100;  % Number of iterations
PFA.lb               = 0; % Lower bound for random feature initialization
PFA.ub               = 1; % upper bound for random feature initialization
PFA.thres            = 0.5; % threshold for feature selection

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end
if ~defaultsfl
    
    if isfield(param,'PFA') 
        if isfield(param.PFA,'N')        && ~isempty( param.PFA.N ),        PFA.c = param.PFA.N; end
        if isfield(param.PFA,'max_Iter') && ~isempty( param.PFA.max_Iter ), PFA.max_Iter = param.PFA.max_Iter; end
        if isfield(param.PFA,'lb')       && ~isempty( param.PFA.lb),        PFA.lb = param.PFA.lb; end
        if isfield(param.PFA,'ub')       && ~isempty( param.PFA.ub ),       PFA.ub = param.PFA.ub; end
        if isfield(param.PFA,'thres')    && ~isempty( param.PFA.thres ),    PFA.thres = param.PFA.thres; end
    end
    nk_PrintLogo
    mestr = 'Path finder algorithm setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
    act = nk_input('Define Path Finder Algorithm (PFA) parameters',0,'mq', ...
                    [sprintf('Number of particles [ %g ]|', PFA.N) ...
                     sprintf('Number of iterations [ %g ]|', PFA.max_Iter) ...
                     sprintf('Lower bound [ %g ]|', PFA.lb) ...
                     sprintf('Upper bound [ %g ]|', PFA.ub) ... 
                     sprintf('Threshold [ %g ]', PFA.thres)], 1:5);
    switch act
        case 1
            PFA.N = nk_input('Define no. of particles',0,'e',PFA.N);
        case 2
            PFA.max_Iter = nk_input('Define maximum no. of iterations',0,'e',PFA.max_Iter);
        case 3
            PFA.lb = nk_input('Define lower bound of feature weight',0,'e',PFA.lb);
        case 4
            PFA.ub = nk_input('Define upper bound of feature weight',0,'e',PFA.ub);
        case 5
            PFA.thres = nk_input('Define threshold for feature selection',0,'e',PFA.thres);
    end
else
    act = 0;
end
param.PFA = PFA; 
