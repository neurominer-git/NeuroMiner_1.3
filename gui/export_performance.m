function [ERR, TBL] = export_performance(handles, batchmode)

warning off 

switch handles.modeflag
    case 'classification'
        sheetstr = 'Classifier'; 
    case 'regression'
        sheetstr = 'Regressor';  
end

filename = sprintf('%s_A%g_PredictionMetrics', ...
        handles.params.TrainParam.SAV.matname, handles.curranal);

for i=1:handles.nclass
    sheetname = sprintf('%s%g',sheetstr,i);
    switch handles.modeflag
        case 'classification'
            TBL = handles.BinClass{i}.tbl_cont;
            if handles.oocvview
                TBL.array=[]; TBL.rownames = [];
                [handles, oocvind] = get_oocvind(handles);
                if isfield(handles.OOCV(oocvind).data,'BinResults')
                    if handles.selSubGroupOOCV.Value>1
                        if iscell(handles.OOCV(oocvind).data.BinResults{handles.curlabel}.Group{handles.selSubGroupOOCV.Value-1}.PredictionPerformance)
                            contigmat = handles.OOCV(oocvind).data.BinResults{handles.curlabel}.Group{handles.selSubGroupOOCV.Value-1}.PredictionPerformance{i};
                        else
                            contigmat = handles.OOCV(oocvind).data.BinResults{handles.curlabel}.Group{handles.selSubGroupOOCV.Value-1}.PredictionPerformance;
                        end
                    else
                        if iscell(handles.OOCV(oocvind).data.BinResults{handles.curlabel}.contingency)
                            contigmat = handles.OOCV(oocvind).data.BinResults{handles.curlabel}.contingency{i};
                        else
                            contigmat = handles.OOCV(oocvind).data.BinResults{handles.curlabel}.contingency;
                        end
                    end
                else
                    contigmat = handles.OOCV(oocvind).data.MultiResults{handles.curlabel}.contingency{i};
                end
                if isfield(handles,'PermAnal'), contigmat.Pvalue = handles.PermAnal; end
            end
        case 'regression'
            TBL = handles.Regr.tbl_cont;
            if handles.oocvview
                [handles, oocvind] = get_oocvind(handles);
                TBL = handles.curRegr.tbl_cont;
                contigmat = handles.curRegr.contigmat;
                if isfield(handles,'PermAnal'), contigmat.Pvalue = handles.PermAnal; end
            end
    end
    if handles.oocvview
        tbl_cont.rownames   = fieldnames(contigmat);
        tbl_cont.array      = struct2cell(contigmat);
        remind = find(strcmp('FPRvec',tbl_cont.rownames) | strcmp('TPRvec', tbl_cont.rownames) | strcmp('X',tbl_cont.rownames));
        tbl_cont.array(remind) = [];
        tbl_cont.array = cell2mat(tbl_cont.array);
        tbl_cont.rownames(remind) = [];
        TBL.colnames = {'Metric', 'Value'};
        TBL.array = [ TBL.array; tbl_cont.array ];
        TBL.rownames = [TBL.rownames; tbl_cont.rownames ];
        filename = sprintf('%s_OOCV%g', filename, oocvind);
        if handles.selSubGroupOOCV.Value>1
            switch handles.modeflag
                case 'classification'
                    groupname = handles.OOCV(oocvind).data.BinResults{handles.curlabel}.Group{handles.selSubGroupOOCV.Value-1}.GroupName;
                case 'regression'
                    groupname = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{handles.selSubGroupOOCV.Value-1}.GroupName;
            end
            filename = sprintf('%s_Group-%s', filename, groupname);
        end
    end
    % Write performance file
    ERR = tbl2file(TBL, filename, sheetname);
end
if isfield(handles,'MultiClass')
    sheetname = sprintf('MultiClass');
    ERR = tbl2file(handles.MultiClass.tbl_cont, filename, sheetname);
end

warning on