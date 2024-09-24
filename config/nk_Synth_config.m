function [SYNTH, act, mess] = nk_Synth_config(SYNTH, mess, parentstr)
global NM

%% Set-up configuration interface

% Define variables
if ~exist('mess','var'), mess = []; end
na_str                      = '?';
synth_flag                  = 2;
methodnum                   = 1;
k                           = 5;
distanceMeasure             = 'Euclidean';
standardize_data            = 2;
numSyntheticObservations    = round(size(NM.label,1)*0.1);
write2disk                  = 1;
adasyn_beta                 = 0.5;
adasyn_kdensity             = 5;
adasyn_kSMOTE               = 5;
numComponents               = 2;
minMethod                   = 1;

% Generate synthetic data?
if isfield(SYNTH,'flag'),                       synth_flag = SYNTH.flag; else, SYNTH.flag = synth_flag; end
if isfield(SYNTH,'method'),                     methodnum = SYNTH.method; else, SYNTH.method = methodnum; end
if isfield(SYNTH,'k'),                          k = SYNTH.k; else, SYNTH.k = k; end
if isfield(SYNTH,'distanceMeasure'),            distanceMeasure = SYNTH.distanceMeasure; else, SYNTH.distanceMeasure = distanceMeasure; end
if isfield(SYNTH,'standardize_data'),           standardize_data = SYNTH.standardize_data; else, SYNTH.standardize_data = standardize_data; end
if isfield(SYNTH,'numSyntheticObservations'),   numSyntheticObservations = SYNTH.numSyntheticObservations; else, SYNTH.numSyntheticObservations = numSyntheticObservations; end
if isfield(SYNTH,'write2disk'),                 write2disk = SYNTH.write2disk; else,  SYNTH.write2disk = write2disk;  end
if isfield(SYNTH,'beta'),                       adasyn_beta = SYNTH.beta; else, SYNTH.beta = adasyn_beta; end
if isfield(SYNTH,'kDensity'),                   adasyn_kdensity = SYNTH.kDensity; else, SYNTH.kDensity = adasyn_kdensity; end            
if isfield(SYNTH,'kSMOTE'),                     adasyn_kSMOTE = SYNTH.kSMOTE; else, SYNTH.kSMOTE = adasyn_kSMOTE; end
if isfield(SYNTH,'n_components'),               numComponents = SYNTH.n_components; else, SYNTH.n_components = numComponents; end
if isfield(SYNTH,'minMethod'),                  minMethod = SYNTH.minMethod; else, SYNTH.minMethod = minMethod ; end

mn_str = []; mn_act = [];

synthstr        = {'enabled','disabled'};
methodstr       = {'kNN-based','ADASYN','PCAGauss','GMM'};
write2diskstr   = {'use data permanently by writing it to disk','recreate data in each model training run'};
stdstr          = {'standardize data', 'data already standardized/standardization not needed'};
mn_str = [mn_str sprintf('Generate synthetic data [ %s ]',synthstr{synth_flag})];       mn_act = [mn_act 1];

if synth_flag == 1
    
    mn_str = [mn_str sprintf('|Choose method for synthetic data generation [ %s ]', methodstr{methodnum})];                                            mn_act = [mn_act 2];
    mn_str = [mn_str sprintf('|Define whether data should be standardized before running distance analysis [ %s ]', stdstr{standardize_data})];        mn_act = [mn_act 3];
    mn_str = [mn_str sprintf('|Choose how to manage synthetic data [ %s ]', write2diskstr{write2disk})];                                               mn_act = [mn_act 4];

    switch methodnum
        case 1
            mn_str = [mn_str sprintf('|Define no. of nearest neighbors for interpolation [ k=%g ]', k)];                                               mn_act = [mn_act 5];
            mn_str = [mn_str sprintf('|Choose distance measure [ %s ]', distanceMeasure)];                                                             mn_act = [mn_act 6];
            mn_str = [mn_str sprintf('|Define how many synthetic observations you would like to be generated [ %g ]', numSyntheticObservations)];      mn_act = [mn_act 7];
            
        case 2
            mn_str = [ mn_str sprintf('|Define beta value (defines how much balancing will be applied, 0<->1) [ beta=%g ]|', adasyn_beta)];             mn_act = [mn_act 8];
            mn_str = [ mn_str sprintf('|Define k value of density algorithm (kNN looking at both classes) [ k=%g ]|',adasyn_kdensity)];                 mn_act = [mn_act 9];
            mn_str = [ mn_str sprintf('|Define k value of SMOTE algorithm (kNN looking only at the minority class) [ k=%g ]|',adasyn_kSMOTE)];          mn_act = [mn_act 10];
            % CV: do you not only need to define numSyntheticObservations
            % for this method? 
        case 3
            mn_str = [mn_str sprintf('|Define how many synthetic observations you would like to generate [ %g ]', numSyntheticObservations)];           mn_act = [mn_act 7];
        case 4
            minMethodStr = {'Bayesian Information Criterion (BIC)', 'Akaike Information Criterion (AIC)'};
            mn_str = [mn_str sprintf('|Define range of gaussian components you want to test in the Gaussian Mixture Modelling [ %s ]', nk_ConcatParamstr(numComponents))]; mn_act = [mn_act 11];
            mn_str = [mn_str sprintf('|Define metric to find optimal number of components [ %s ]', minMethodStr{minMethod})];                           mn_act = [mn_act 12];
            mn_str = [mn_str sprintf('|Define how many synthetic observations you would like to generate [ %g ]', numSyntheticObservations)];           mn_act = [mn_act 7];

    end
end

nk_PrintLogo

if ~isempty(mess)
    for i=1:numel(mess)
        if isempty(mess(i).text), continue; end
        fprintf('\n');mess(i).text = regexprep(mess(i).text,'\','/');
        fprintf(mess(i).text); 
    end
    fprintf('\n')
    mess = [];
end

fprintf('\n'); mestr = 'Configure synthetic data generation';  
navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); fprintf('You are here: %s >>> ',parentstr); 
act = char(nk_input(mestr,0,'mq', mn_str, mn_act));

switch act
    
    case 1
       if SYNTH.flag == 1, SYNTH.flag = 2; else, SYNTH.flag = 1; end
    case 2
        SYNTH.method = nk_input('Choose synthetic data generation method', 0, 'm', 'kNN-based|ADASYN|PCAGauss|GMM', 1:4, methodnum);
    case 3
        if SYNTH.standardize_data == 1, SYNTH.standardize_data = 2; else, SYNTH.standardize_data = 1; end

    case 4
        if SYNTH.write2disk == 1, SYNTH.write2disk = 2; else, SYNTH.write2disk = 1; end
    
    case 5
        SYNTH.k = nk_input('Define no. of nearest neighbors for interpolation',0,'i', k);
    
    case 6
       sel = 'Euclidean|Manhattan|Cosine|Mahalanobis';
       defsel = strsplit(sel,'|');
       SYNTH.distanceMeasure = char(nk_input('Choose distance measure',0,'m', sel, defsel, find(contains(defsel, distanceMeasure))));
    
    case 7
        SYNTH.numSyntheticObservations = nk_input('How many synthetic observations should be generated',0,'i', numSyntheticObservations);

    case 8
        SYNTH.beta = nk_input('Define beta for degree of balancing',0,'e',adasyn_beta);
        
    case 9
        SYNTH.kDensity = nk_input('Define k for density estimation in ADASYN',0,'i',adasyn_kDensity);
        
    case 10
        SYNTH.kSMOTE = nk_input('Define k for SMOTE in ADASYN',0,'i',adasyn_kSMOTE);

    case 11
        SYNTH.n_components = nk_input('Define range of GMM component numbers',0,'i',numComponents);
    
    case 12
        if SYNTH.minMethod ==1, SYNTH.minMethod = 2; else, SYNTH.minMethod = 1; end  

end