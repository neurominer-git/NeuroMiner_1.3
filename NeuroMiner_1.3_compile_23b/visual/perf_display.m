function handles = perf_display(handles)
% Set current analysis
analind = get(handles.selAnalysis,'Value'); 
if ~isfield(handles,'modeflag')
    % set alternative label
    if isfield(handles.NM.analysis{1,analind}.params,'label')
        handles.label = handles.NM.analysis{1,analind}.params.label.label;
        handles.modeflag = handles.NM.analysis{1,analind}.params.label.modeflag;
    else
        handles.label = handles.NM.label;
        handles.modeflag = handles.NM.modeflag;
    end
end
if ~isfield(handles,'curranal')
    handles.curranal = 1; 
end

handles.prevanal = handles.curranal;
handles.curranal = analind;

if ~handles.NM.analysis{analind}.status || ~isfield(handles.NM.analysis{analind},'GDdims')
    set_visibility(handles)
    set_panel_visibility(handles,'off');
    return; 
end

% Set current modality
% Set multi-modal flag
if numel(handles.NM.analysis{analind}.GDdims)>1
    handles.multi_modal = 1;
    set(handles.selModal,'Enable','on')
    load_modal(handles, handles.NM.analysis{analind}.GDdims);
    if isfield(handles.NM.analysis{analind},'META')
       popuplist = handles.selModal.String;
       popuplist{end+1}='Bagged predictor';
       handles.selModal.String = popuplist;
    end
else
    handles.multi_modal = 0;
    set(handles.selModal,'Enable','off')  
end

if isfield(handles,'multilabel') && handles.multilabel
    if isfield(handles.NM.analysis{analind}.params.TrainParam,'MULTILABEL')
        if isfield(handles,'curlabel') && handles.curlabel > numel(handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel)
            handles.selLabel.Value = 1;
        end
        handles.selLabel.String = handles.NM.labelnames(handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel);
        handles.curlabel = handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel(get(handles.selLabel,'Value'));
    else
        handles.selLabel.String = handles.NM.labelnames;
        handles.curlabel= get(handles.selLabel,'Value');
    end
end

[handles, visdata] = switch_analysis(handles);

if isfield(handles,'MLIapp') && ~isfield(handles.NM.analysis{analind},'MLI') && ~isempty(handles.MLIapp) && handles.MLIapp~=0
    handles.MLIapp.delete;
    handles = rmfield(handles,'MLIapp');
elseif isfield(handles,'MLIapp') && handles.MLIapp~=0 && isfield(handles.NM.analysis{analind},'MLI') && ishandle(handles.MLIapp)
    updateFcn(handles.MLIapp,handles);
end

handles.lbStartup.String = 'Customize menus ...';
set_panel_visibility(handles,'on')
set_visibility(handles)
load_selYAxis(handles)
load_popupmenu1(handles)
load_selModelMeasures(handles)
load_selSubParams(handles)

if ~isempty(visdata)
    load_selModality(handles); 
    load_selPager(handles); 
else
    handles.selModelMeasures.Value = 1;
end

if isfield(handles,'MultiClass'), load_selOneVsAll_Info(handles); end
handles = display_main(handles);





