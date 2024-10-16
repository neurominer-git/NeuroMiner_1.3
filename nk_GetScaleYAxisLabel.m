function [yl, ylb, yfun, ylbshort] = nk_GetScaleYAxisLabel(param)
global PREPROC NM MULTILABEL

if ~exist('param','var')
    OPTCRIT = NM.TrainParam.SVM.GridParam;
elseif isnumeric(param)
    OPTCRIT = param;
else
    OPTCRIT = param.GridParam;
end
    
switch OPTCRIT
    case 1
        yl = [0 100];   ylb = 'Accuracy [%]';                               yfun = @ACCURACY;    ylbshort = 'ACCURACY';
    case 2
        yl = [0 100];   ylb = 'Sensitivity [%]';                            yfun = @SENSITIVITY; ylbshort = 'SENSITIVITY';
    case 3
        yl = [0 100];   ylb = 'Specificity [%]';                            yfun = @SPECIFICITY; ylbshort = 'SPECIFICITY';
    case 4
        yl = [0 100];   ylb = 'False Positive Rate [%]';                    yfun = @FPR; ylbshort = 'FPR';
    case 5
        yl = [0 100];   ylb = 'Positive Predictive Value [%]';              yfun = @PPV; ylbshort = 'PPV';
    case 6
        yl = [-1 1];    ylb = 'Matthews Correlation Coefficient';           yfun = @MCC; ylbshort = 'MCC';      
    case 7
        yl = [0 1];     ylb = 'Area-Under-the-Curve';                       yfun = @AUC; ylbshort = 'AUC';
    case 9
        yl = [0 range(NM.label(:,MULTILABEL.dim)).^2];   ylb = 'Mean squared error'; yfun = @MSE; ylbshort = 'MSE';
    case 10
        yl = [0 100];   ylb = 'Squared correlation coefficient [% explained variance]'; yfun = @SCC; ylbshort = 'SCC';
    case 11
        yl = [0 100];   ylb = 'Normalized root of mean squared deviation [%]'; yfun = @NRMSD; ylbshort = 'NRMSD';
    case 12
        yl = [0 1];     ylb = 'Root of mean squared deviation';             yfun = @RMSD; ylbshort = 'RMSD';
    case 13
        yl = [0 1];     ylb = 'Gmean';                                      yfun = @GMEAN; ylbshort = 'GMEAN';
    case 14
        yl = [0 100];   ylb = 'Balanced Accuracy [%]';                      yfun = @BAC; ylbshort = 'BAC';
    case 15
        yl = [0 100];   ylb = 'F1-Score';                                   yfun = @FSCORE; ylbshort = 'FSCORE';
    case 16
        yl = [-1 1];    ylb = 'correlation coefficient';                    yfun = @CC; ylbshort = 'CC';
    case 17
        yl = [0 100];   ylb = 'Enhanced Balanced Accuracy [%]';             yfun = @BAC2; ylbshort = 'BAC2';
    case 19
         yl = [-100 100];   ylb = 'Prognostic Summary Index [%]';           yfun = @PSI; ylbshort = 'PSI';
    case 20
         yl = [0 10]; ylb = 'Number Needed to Predict';                     yfun = @NNP; ylbshort = 'NNP';
    case 21
         yl = [0 1]; ylb = 'Expected Calibration Error';                    yfun = @ECE; ylbshort = 'ECE'; 
    case 22
         yl = [0 15]; ylb = 'Positive Likelihood Ratio';                  yfun = @pLR; ylbshort = 'pLR'; 
    case 23
         yl = [0 1]; ylb = 'Negative Likelihood Ratio';                  yfun = @nLR; ylbshort = 'nLR'; 
    case 18
        if iscell(PREPROC)
            iPREPROC = PREPROC{1}; else, iPREPROC = PREPROC;
        end
        if isfield(iPREPROC,'LABELMOD') && ...
            isfield(iPREPROC.LABELMOD,'TARGETSCALE') && ...
            iPREPROC.LABELMOD.TARGETSCALE
            yl = [0 1]; 
        else
            if ~isempty(NM)
                yl = [0 nk_Range(NM.label) ];
            else
                yl = [0 nk_Range(evalin('base','NM.label'))];
            end
        end
        ylb = 'Mean Average Error [Label Range]';                           yfun = @MAERR; ylbshort = 'MAE';
end

