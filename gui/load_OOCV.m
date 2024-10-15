function handles = load_OOCV(handles, oocv)

cnt=1; l = handles.selLabel.Value;

% Loop through all OOCV containers
for n = 1: numel(oocv)
    
   % If analysis has been applied to OOCV container, retrieve predictions
   % and prepare export tables 
   if ~isempty(oocv{n})

        handles.OOCV(cnt).num = n;
        handles.OOCV(cnt).data = oocv{n};
        
        % Do we have to deal with subgroups?
        if isfield(handles.NM.OOCV{n},'groups')
            handles.OOCV(cnt).groups = handles.NM.OOCV{n}.groups;
            handles.OOCV(cnt).grpnames = handles.NM.OOCV{n}.grpnames;
            if isfield(handles.NM.OOCV{n},'refgroup')
                handles.OOCV(cnt).refgroup = handles.NM.OOCV{n}.refgroup;
            end
        end

        % Do we know the OOCV data's labels
        if handles.OOCVinfo.Analyses{handles.curranal}.labels_known(cnt)
            labels_known = true; 
        else 
            labels_known = false; 
        end

        % Create table data for export
        % Handle the Regression scenario differently
        if isfield(handles,'Regr')
            handles.OOCV(cnt).data.tbl = struct('rownames',[],'colnames',[],'array',[]);
            handles.OOCV(cnt).data.tbl.rownames = handles.OOCVinfo.Analyses{handles.curranal}.cases{n};  
            % Are the labels known or not?
            if labels_known
                handles.OOCV(cnt).data.tbl.colnames = {'Cases', 'EXP_LABEL', 'Mean_Pred', 'Std_Pred'};
                handles.OOCV(cnt).data.tbl.array    = handles.OOCVinfo.Analyses{handles.curranal}.label{n};
            else
                handles.OOCV(cnt).data.tbl.colnames = {'Cases', 'Mean_Pred', 'Std_Pred'};
                handles.OOCV(cnt).data.tbl.array    = [];
            end
            handles.OOCV(cnt).data.tbl.array = [ handles.OOCV(cnt).data.tbl.array ...
                                                 handles.OOCV(cnt).data.RegrResults{l}.MeanCV2PredictedValues ...
                                                 handles.OOCV(cnt).data.RegrResults{l}.StdCV2PredictedValues ] ;
            if isfield(handles.OOCV(cnt).data,'RegrResults')
                handles.OOCV(cnt).flag = true;
            else
                handles.OOCV(cnt).flag = false;
                cnt=cnt+1;
                continue
            end
        else
            nclass = numel(handles.BinClass);
            handles.OOCV(cnt).data.tbl = struct('rownames',[],'colnames',[],'array',[]);
            if isfield( handles.OOCV(cnt).data,'MultiResults'), multi_res = true; else, multi_res=false; end
            if multi_res
                Cont = handles.OOCV(cnt).data.MultiResults{l};
                handles.OOCV(cnt).flag = true;
            else
                if isfield(handles.OOCV(cnt).data,'BinResults')
                    if iscell(handles.OOCV(cnt).data.BinResults)
                        Cont = handles.OOCV(cnt).data.BinResults{1};
                    else
                        Cont = handles.OOCV(cnt).data.BinResults;
                    end
                    handles.OOCV(cnt).flag = true;
                else
                    handles.OOCV(cnt).flag = false;
                    cnt=cnt+1;
                    continue
                end
            end
            if handles.OOCV(cnt).flag
                for j=1:nclass
                    % Cases are always the same in the OOCV setting
                    handles.OOCV(cnt).data.tbl(j).rownames = handles.OOCVinfo.Analyses{handles.curranal}.cases{cnt};    
                    if labels_known
                        handles.OOCV(cnt).data.tbl(j).colnames = {'Cases', 'EXP_LABEL', 'PRED_LABEL', 'Errors', 'Mean_Score', 'Std_Score', 'Ens_Prob1', 'Ens_Prob-1', 'PercRank_Score'};
                        handles.OOCV(cnt).data.tbl(j).array    = Cont.BinLabels{j}; Errors =  Cont.BinProbPredictions{j} ~=  Cont.BinLabels{j};
                    else
                        handles.OOCV(cnt).data.tbl(j).colnames = {'Cases', 'PRED_LABEL', 'Mean_Score', 'Std_Score', 'Ens_Prob1', 'Ens_Prob-1', 'PercRank_Score'};
                        handles.OOCV(cnt).data.tbl(j).array    = []; Errors = [];
                    end
                    
                    handles.OOCV(cnt).data.tbl(j).array = [ handles.OOCV(cnt).data.tbl(j).array ...
                                                             Cont.BinProbPredictions{j} ...
                                                             Errors ...
                                                             Cont.MeanCV2PredictedValues{j} ...
                                                             Cont.StdCV2PredictedValues{j} ...
                                                             Cont.BinMajVoteProbabilities{j} ...
                                                             ranktransform(handles.BinClass{j}.mean_predictions, Cont.MeanCV2PredictedValues{j})] ; 
                end
                if multi_res
                    handles.OOCV(cnt).data.tbl_mult = struct('rownames',[],'colnames',[],'array',[]);
                    handles.OOCV(cnt).data.tbl_mult.rownames = handles.OOCVinfo.Analyses{handles.curranal}.cases{cnt};    
                    if labels_known
                        handles.OOCV(cnt).data.tbl_mult.colnames = {'Cases', 'EXP_LABEL', 'PRED_LABEL', 'Errors'};
                        handles.OOCV(cnt).data.tbl_mult.array    = handles.OOCVinfo.Analyses{handles.curranal}.label{cnt}; Errors = Cont.MultiClass.multi_predictionsCV2 ~= handles.OOCVinfo.Analyses{handles.curranal}.label{cnt};
                    else
                        handles.OOCV(cnt).data.tbl_mult.colnames = {'Cases', 'PRED_LABEL' };
                        handles.OOCV(cnt).data.tbl_mult.array    = []; Errors = [];
                    end
                    for j=1:nclass
                         handles.OOCV(cnt).data.tbl_mult.colnames{end+1} = sprintf('Prob-G%g',j);
                    end
                    handles.OOCV(cnt).data.tbl_mult.array = [ handles.OOCV(cnt).data.tbl_mult.array ...
                                                             Cont.MultiClass.multi_predictionsCV2 ...
                                                             Errors ...
                                                             Cont.MultiClass.multi_probabilitiesCV2] ;
                end
            end
        end
        cnt = cnt + 1;
   end
end