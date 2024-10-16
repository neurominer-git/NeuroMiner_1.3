function handles = display_main(handles)

h_class         = handles.popupmenu1.Value;
h_classlist     = handles.popupmenu1.String;
h_list          = handles.selModelMeasures.String;
h_val           = handles.selModelMeasures.Value;
if h_val>numel(h_list), h_val=1; end
handles.curclass = h_class; 

if isfield(handles,'SubIndex') 
    handles = rmfield(handles,'SubIndex');
end

switch h_list{h_val}
    
    case 'Classification plot'
     
        handles.pnModelPerf.Visible         = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnBinary.Visible            = 'on';
  
        if strcmpi(h_classlist{h_class},'Multi-group classifier')
            handles.oocvview = false;
            handles.cmdExportCobWeb.Visible = 'on';
            handles.selOneVsAll_Info.Enable  = 'on';
            %set(get(handles.pnBinRegrPerfCmd, 'Children'), 'Enable', 'off');
            load_selYAxis(handles)
            load_selModelMeasures(handles)
            
            handles = display_multiclassplot(handles);
            handles = sel_onevsone(handles, handles.selOneVsAll_Info);
            load_selCase(handles,handles.MultiClass.cases)
        else
            %set(get(handles.pnBinRegrPerfCmd, 'Children'), 'Enable', 'on');
            if isfield(handles,'MultiClass') && isfield(handles.MultiClass,'spideraxes') 
                handles.MultiClass.spideraxes.Title.Visible='off';
            end
            handles.selOneVsAll_Info.Enable = 'off';
            handles.cmdExportCobWeb.Visible = 'off';
            handles.cmdMetricExport.Visible = 'on';
            load_selYAxis(handles)
            load_selModelMeasures(handles)
            
            if strcmp(handles.selCVoocv.Enable,'on') && handles.selCVoocv.Value>1
                handles.oocvview = true;
                [~,oocvind] = get_oocvind(handles);
                if handles.OOCV(oocvind).flag
                    if isfield(handles,'MultiClass'), fldname = 'MultiResults'; else, fldname = 'BinResults'; end
                    if isfield(handles.OOCV(oocvind).data.(fldname){handles.curlabel},'PermAnal')
                        handles.PermAnal = handles.OOCV(oocvind).data.(fldname){handles.curlabel}.PermAnal.ModelPermSignificance(h_class);
                    else
                        if isfield(handles,'PermAnal'), handles = rmfield(handles,'PermAnal'); end
                    end
                    handles = display_classplot_oocv(h_class, handles);
                    load_selCase(handles,handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind});
                    if isfield(handles.OOCV(oocvind).data,fldname) && isfield(handles.OOCV(oocvind).data.(fldname){h_class},'Group')
                        Groups = handles.OOCV(oocvind).data.(fldname){h_class}.Group;
                        GroupNames = cell(numel(Groups)+1,1);
                        GroupNames{1} = 'Show entire OOCV sample';
                        for g=2:numel(Groups)+1
                            GroupNames{g} = Groups{g-1}.GroupName;
                        end
                        handles.selSubGroupOOCV.String = GroupNames;
                        handles.selSubGroupOOCV.Value = 1;
                        handles.selSubGroupOOCV.Visible = 'on';                    
                    else
                        handles.selSubGroupOOCV.Value = 1;
                        handles.selSubGroupOOCV.Visible = 'off';
                    end
                end
            else
                handles.oocvview = false;
                if isfield(handles,'PermAnal'), handles = rmfield(handles,'PermAnal'); end
                handles = display_classplot(h_class, handles);
                load_selCase(handles,handles.BinClass{h_class}.cases)
                handles.selSubGroupOOCV.Value = 1;
                handles.selSubGroupOOCV.Visible = 'off';
            end
        end
        
    case 'Regression plot'
        
        handles.pnModelPerf.Visible         = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnBinary.Visible            = 'on';
        handles.cmdExportCobWeb.Visible     = 'off';
        handles.cmdMetricExport.Visible     = 'on';
        cla(handles.axes20); 
        handles.axes20.Visible='off'; drawnow
        handles.cmdExportAxes20.Visible = 'off';
        

        if strcmp(handles.selCVoocv.Enable,'on') && handles.selCVoocv.Value > 1

            [~,oocvind] = get_oocvind(handles);
            if isfield(handles.OOCV(oocvind).data.RegrResults{1},'Group')
                Groups = handles.OOCV(oocvind).data.RegrResults{1}.Group;
                GroupNames = cell(numel(Groups)+1,1);
                GroupNames{1} = 'Show entire OOCV sample';
                for g=2:numel(Groups)+1
                    GroupNames{g} = Groups{g-1}.GroupName;
                end
                handles.selSubGroupOOCV.String = GroupNames;
                handles.selSubGroupOOCV.Value = 1;
                handles.selSubGroupOOCV.Visible = 'on';
            else
                handles.selSubGroupOOCV.Value = 1;
                handles.selSubGroupOOCV.Visible = 'off';
            end

            [handles, oocvind] = get_oocvind(handles);
            % Check whether the labels are known
            labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);

            if isfield(handles.OOCV(oocvind).data.RegrResults{handles.curlabel},'PermAnal')
                handles.PermAnal = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.PermAnal.ModelPermSignificance;
            else
                if isfield(handles,'PermAnal'), handles = rmfield(handles,'PermAnal'); end
            end
            
            if ~labels_known
                handles = display_regrplot_oocv_labels_unknown(h_class, handles);
            else
                handles  = display_regrplot(handles, [], false, true, false, 0.2);
                handles.oocvview = true;
                handles  = display_regrplot(handles, [], handles.oocvview, false, true, 0.8);
                load_selCase(handles,handles.OOCVinfo.Analyses{handles.curranal}.cases{handles.oocvind});
            end
          
        else
            handles.oocvview = false;
            if isfield(handles,'PermAnal'), handles = rmfield(handles,'PermAnal'); end
            handles = display_regrplot(handles, [], false, false, true, 0.8);
            handles.selSubGroupOOCV.Value = 1;
            handles.selSubGroupOOCV.Visible = 'off';
            load_selCase(handles,handles.Regr.cases)
        end

    case 'Visualization results'
      
        handles.pnModelPerf.Visible         ='off';
        handles.pnBinary.Visible            ='off';
        handles.pnVisual.Visible            ='on';
        load_selModality(handles)
        handles = display_visual(handles);

    case 'ML Interpreter results'
        
        %handles.pnModelPerf.Visible         ='off';
        %handles.pnBinary.Visible            ='off';
        %handles.pnVisual.Visible            ='on';
        %load_selModality(handles)
        handles.MLIapp = appMLI(handles);
        handles.selModelMeasures.Value = 1;
        display_main(handles);

    otherwise
      
        handles.pnBinary.Visible            = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnModelPerf.Visible         = 'on';
        load_selModelMeasures(handles)
        switch handles.METAstr
            case 'none'
                handles                     = display_modelperf(handles);
            case 'bagged'
                handles                     = display_modelperf_bagged(handles);
        end
end
