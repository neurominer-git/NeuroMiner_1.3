function [act, param] = nk_LIBLIN_config(modeflag, param, defaultsfl, parentstr)
% Nikolaos Koutsouleris: Added bias term configuration option on the 19/02/2023
classifier = 0;
tolerance = 0.01;
b0 = 0;
weighting = 0;
bias = -1;
AutoFl = false; autoflag = false; autostr = 'no'; 
NonRegBias = false; NonRegBiasFlag = false; NonRegBiasStr = 'no'; 

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end

% Type of classifier
if ~exist('param','var') || ~isfield(param,'LIBLIN') 
    param.LIBLIN = []; 
else
    if isfield(param.LIBLIN,'classifier')
        classifier = param.LIBLIN.classifier;
    end
    if isfield(param.LIBLIN,'tolerance')
        tolerance = param.LIBLIN.tolerance;
    end
end

% Weighting of the hyperplane
if isfield(param.LIBLIN,'Weighting')
    weighting = param.LIBLIN.Weighting;
end
switch weighting
    case 1
        weightstr = 'yes';
    case 0
        weightstr = 'no';
end

% Bias definition
if isfield(param.LIBLIN,'bias')
    bias = param.LIBLIN.bias;
end
switch bias
    case 1
        biasstr = 'yes';
    case -1
        biasstr = 'no';
end

b0flag = false;

switch modeflag
    case 'classification'
        switch classifier
            case 0
                cldef = 1; clstr = 'L2-regularized logistic regression (primal)';
                b0flag = true;
                autoflag = true;
                NonRegBiasFlag = true;
            case 1
                cldef = 2; clstr = 'L2-regularized L2-loss support vector classification (dual)';
            case 2
                cldef = 3; clstr = 'L2-regularized L2-loss support vector classification (primal)';
                autoflag = true;
                NonRegBiasFlag = true;
            case 3
                cldef = 4; clstr = 'L2-regularized L1-loss support vector classification (dual)';
            case 4

            case 5
                cldef = 5; clstr = 'L1-regularized L2-loss support vector classification';
                NonRegBiasFlag = true;
            case 6
                cldef = 6; clstr = 'L1-regularized logistic regression';
                b0flag = true;
                NonRegBiasFlag = true;
            case 7
                cldef = 7; clstr = 'L2-regularized logistic regression (dual)';
                b0flag = true;
            otherwise
                cldef=1; clstr = 'NA';
        end
    case 'regression'
        switch classifier
            case 11
                cldef=1; clstr = 'L2-regularized L2-loss epsilon support vector regression (primal)';
                autoflag = true;
                NonRegBiasFlag = true;
            case 12 
                cldef=2; clstr = 'L2-regularized L2-loss epsilon support vector regression (dual)';
            case 13
                cldef=3; clstr = 'L2-regularized L1-loss epsilon support vector regression (dual)';
            otherwise 
                cldef=1; clstr = 'NA';
        end
                
end

% Probability estimates
if isfield(param.LIBLIN,'b')
    switch param.LIBLIN.b
        case 0
            b0str = 'no'; 
        otherwise
            b0str = 'yes';
    end
    b0 = param.LIBLIN.b;
else
    param.LIBLIN.b = 0; b0str = 'no'; 
end

if isfield(param.LIBLIN,'AutoC') && autoflag
    AutoFl = param.LIBLIN.AutoC;
    switch AutoFl
        case 0
            autostr = 'no'; 
        otherwise
            autostr = 'yes';
    end
else
    AutoFl = false;
end

if isfield(param.LIBLIN,'NonRegBias') && NonRegBiasFlag && bias == 1
    NonRegBias = param.LIBLIN.NonRegBias;
    switch NonRegBias
        case 0
            NonRegBiasStr = 'no'; 
        otherwise
            NonRegBiasStr = 'yes';
    end
else
    NonRegBias = false;
end

if ~defaultsfl
    switch modeflag
        case 'classification'
            if b0flag
                b0opt = ['Output probability scores [ ' b0str ' ]|']; b0num = 3;
            else
                b0opt = ''; b0num = [];
            end
            menuact = [ 'Classifier type [ ' clstr ' ]|' ...
                        'Weighting of hyperplane in unbalanced problems [ ' weightstr ' ]|' ...
                        b0opt ...
                        'Tolerance [ ' num2str(tolerance) ' ]|'];
            menusel = [1 2 b0num 4];
           
        case 'regression'
            menuact = [ 'Regressor type [ ' clstr ' ]|' ...
                        'Tolerance [ ' num2str(tolerance) ' ]|'];
            menusel = [1 4];
                       
    end
    menuact = [ menuact 'Use bias term [ ' biasstr ' ]' ];
    menusel = [ menusel 5 ];
    if NonRegBiasFlag && bias == 1
        menuact = [ menuact '|Do not regularize bias [ ' NonRegBiasStr ' ]' ];
        menusel = [ menusel 6 ];
    end
    if autoflag
        menuact = [ menuact '|Auto-optimize LIBLINEAR hyperparameters (make sure you set them to 1) [ ' autostr ' ]' ];
        menusel = [ menusel 7 ];
    end

    nk_PrintLogo
    mestr = 'LIBLINEAR: Parameter setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
    act = nk_input(mestr ,0 , 'mq', menuact, menusel);
    switch act    
        case 1
            switch modeflag 
                case 'classification'
                    classifier = nk_input('Type of classifier: ',0, 'mq', ...
                                ['L2-regularized logistic regression (primal)|' ...
                                 'L2-regularized L2-loss SVC (dual, similar to LIBSVM)|' ...
                                 'L2-regularized L2-loss SVC (primal, faster than dual implementation)|' ...
                                 'L2-regularized L1-loss SVC (dual, less sensitive to outliers than L2-loss SVC)|' ...
                                 'L1-regularized L2-loss SVC (sparse SVC algorithm for embedded feature selection)|' ...
                                 'L1-regularized logistic regression (sparse LR algorithm for embedded feature selection)|' ...
                                 'L2-regularized logistic regression (dual)'], ...
                                 [1, 2, 3, 4, 6, 7, 8], cldef);
                case 'regression'
                    classifier = nk_input('Type of predictor: ',0, 'mq', ...
                            ['L2-regularized L2-loss SVR (primal, faster than dual implementation)|' ...
                             'L2-regularized L2-loss SVR (dual)|' ...
                             'L2-regularized L1-loss SVR (dual, less sensitive to outliers than L2-loss SVR)|'], ...
                             [12, 13, 14], cldef);
                    weighting = 0;
            end
            if classifier ~= 7 || classifier ~= 1, b0=0; end
            if classifier, classifier = classifier-1; end
        case 2   
            if weighting , weighting = 0; else, weighting = 1; end
        case 3
            if b0, b0 = 0; else, b0 = 1; end
        case 4
            tolerance = nk_input('Set tolerance of termination criterion',0,'e',tolerance);
        case 5
            if bias == 1, bias = -1; else, bias =1; end
        case 6
            NonRegBias = ~NonRegBias; 
        case 7
            AutoFl = ~AutoFl;
    end
end
param.LIBLIN.classifier = classifier;
param.LIBLIN.Weighting  = weighting;
param.LIBLIN.b          = b0;
param.LIBLIN.tolerance  = tolerance;
param.LIBLIN.quiet      = 1;
param.LIBLIN.bias       = bias;
param.LIBLIN.NonRegBias = NonRegBias;
param.LIBLIN.AutoC      = AutoFl;
end