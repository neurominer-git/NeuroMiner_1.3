% =========================================================================
% =                        CLASSIFICATION PLOTS                           =
% =========================================================================
function handles = display_classplot_oocv(h, handles)

% Prepare axes
warning off
axes(handles.axes1); cla(handles.axes1); hold on
handles.axes1.Position = handles.axes1pos_orig;
handles.axes38.Visible = 'off';  cla(handles.axes38);
handles.tglSort.Enable =  "off";
MS = 15;
MSoocv = 40;

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

[handles, oocvind] = get_oocvind(handles);

% Check whether the labels are known
labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);

clrswp = handles.tglClrSwp.Value;
if clrswp
    g1 = 2; g2 = 1;
else
    g1 = 1; g2 = 2;
end

% Current label
l=1;

% Determine if we have a multi-class analysis
if isfield(handles,'MultiClass')
    P_fld = 'MultiResults';
else
    P_fld = 'BinResults';
end

if handles.tglPercRank.Value
    switch GraphType
        case {1,2,3}
            refsample =  handles.BinClass{h}.mean_predictions;
            targetsample = handles.OOCV(oocvind).data.(P_fld){l}.MeanCV2PredictedValues{h};
            zeroline = 0;
        case 4
            refsample =  handles.BinClass{h}.prob_predictions(:,1);
            targetsample = handles.OOCV(oocvind).data.(P_fld){l}.BinMajVoteProbabilities{h};
            zeroline = 0;
        case {5,6}
            refsample =  handles.BinClass{h}.mean_predictions;
            targetsample = handles.OOCV(oocvind).data.(P_fld){l}.BinMajVoteProbabilities{h};
            zeroline = 0;
    end
    zeroln = nk_ComputePercentiles(refsample, zeroline,'inverse');
    P_h = ranktransform(refsample);
    P_oocv_h = ranktransform(refsample, targetsample);
else
    switch GraphType
    case {1,2,3}
        P_h         = handles.BinClass{h}.mean_predictions;
        P_oocv_h    = handles.OOCV(oocvind).data.(P_fld){l}.MeanCV2PredictedValues{h};
        zeroln = 0;
    case 4
        % Majority voting probabilities
        P_h = handles.BinClass{h}.prob_predictions(:,1);
        P_oocv_h    = handles.OOCV(oocvind).data.(P_fld){l}.BinMajVoteProbabilities{h};
        zeroln = 0.5;
    case {5,6}
        % Mean Majority voting probabilities
        P_h = handles.BinClass{h}.CV2grid.mean_predictions;
        P_oocv_h    = handles.OOCV(oocvind).data.(P_fld){l}.BinMajVoteProbabilities{h};
        zeroln = 0;
    end
    if handles.tglP.Value
        P_oocv_h = probtransform(handles.BinClass{h}.labelh,P_h,P_oocv_h,'platt');
        P_h = probtransform(handles.BinClass{h}.labelh,P_h,[],'platt');
        zeroln = 0.5;
    end
end

% Get subindex if availableSubIndex
if isfield(handles,'SubIndex'), subfl = true; SubI = handles.SubIndex; else, subfl = false; SubI = true(size(P_oocv_h,1),1); end

% Mark groups with color
if isfield(handles.BinClass{h},'ind2')
    id1 = handles.BinClass{h}.ind1; 
    id2 = handles.BinClass{h}.ind2; 
else
    id1 = handles.BinClass{h}.ind1; 
    id2 = ~handles.BinClass{h}.ind1; 
end

legvecn = false(1,5);
handlevecn = cell(1,5);

% Print CV data: Group 1
violin(P_h(id1),'x',1,'facecolor',handles.colptin(handles.BinClass{h}.groupind(g1),:),'edgecolor','none', 'facealpha', 0.1, 'bw',0.3);% 'mc','','medc','');
handlevecn{1} = dotdensity( 1 ,P_h(id1), ...
    'dotEdgeColor', handles.colptin(handles.BinClass{h}.groupind(g1),:), ...
    'dotFaceColor',handles.colptin(handles.BinClass{h}.groupind(g1),:), ...
    'dotSize',MS, ...
    'dotMarker','o', ...
    'dotAlpha', 0.5, ...
    'meanLine', 'off', ...
    'medianLine', 'off'); 

legvecn(1)=1;

% Print CV data: Group 2
if numel(handles.BinClass{h}.groupind) == 2
    if ~handles.BinClass{h}.one_vs_all
        CLP = 'o'; CLR = handles.colptin(handles.BinClass{h}.groupind(g2),:);
        handles.tglClrSwp.Enable = "on";
    else
        CLP = 'o'; CLR = rgb('DarkGrey');
        handles.tglClrSwp.Enable = "off";
    end
else
    CLP = 'o'; CLR = 'k';
end

violin(P_h(id2),'x',3,'facecolor',handles.colptin(handles.BinClass{h}.groupind(g2),:),'edgecolor','none', 'facealpha', 0.1, 'bw',0.3);% 'mc','','medc','');
handlevecn{2} = dotdensity(3,P_h(id2), ...
    'dotEdgeColor', CLR, ...
    'dotFaceColor', CLR, ...
    'dotSize',MS, ...
    'dotMarker', CLP, ...
    'dotAlpha', 0.5, ...
    'meanLine', 'off', ...
    'medianLine', 'off');

legvecn(2)=1;

if labels_known
    label_oocv_h = handles.OOCV(oocvind).data.(P_fld){l}.BinLabels{h};
    ind0 =  label_oocv_h ~=0 ; fid0_oocv = find(~label_oocv_h);
    tP_oocv_h = P_oocv_h(ind0 & SubI);
    id1_oocv = label_oocv_h(ind0 & SubI) == 1 & ~isnan(tP_oocv_h) ; fid1_oocv = find(label_oocv_h(ind0 & SubI) == 1);
    id2_oocv = label_oocv_h(ind0 & SubI) == -1 & ~isnan(tP_oocv_h) ; fid2_oocv = find(label_oocv_h(ind0 & SubI) == -1);
    if sum(id1_oocv), oocvfl1 = true; else, oocvfl1 = false; end
    if sum(id2_oocv), oocvfl2 = true; else, oocvfl2 = false; end
    kpos1 = 2; kpos2 = 2; 
    if (oocvfl1 && oocvfl2) 
        kpos1 = 2; kpos2 = 4;
    end
     % Print independent sample prediction using dot density plots: Group 1
     if sum(id1_oocv)
        violin(tP_oocv_h(id1_oocv),'x',kpos1 ,'facecolor',handles.colptin(handles.BinClass{h}.groupind(g1),:),'edgecolor','none', 'facealpha', 0.3, 'bw',0.3);% 'mc','','medc','');
        [handlevecn{3},~,N{1},X{1},Y{1}] = dotdensity( kpos1 ,tP_oocv_h(id1_oocv), ...
            'dotEdgeColor', handles.colptin(handles.BinClass{h}.groupind(g1),:), ...
            'dotFaceColor', handles.colptin(handles.BinClass{h}.groupind(g1),:), ...
            'dotSize',MSoocv, ...
            'dotMarker','o', ...
            'meanLine', 'off', ...
            'medianLine', 'off');
        legvecn(3) = 1;
        fN{1} = fid1_oocv(N{1});
     end
    
     % Print independent sample prediction using dot density plots: Group 2
    if sum(id2_oocv)
        if numel(handles.BinClass{h}.groupind) == 2
            if ~handles.BinClass{h}.one_vs_all
                CLP = 'o'; CLR = handles.colptin(handles.BinClass{h}.groupind(g2),:);
            else
                 CLP = 'o'; CLR = rgb('DarkGrey');
            end
        else
            CLP = 'o'; CLR = 'k';
        end
        violin(tP_oocv_h(id2_oocv),'x',kpos2 ,'facecolor',CLR,'edgecolor','none', 'facealpha', 0.3, 'bw',0.3);% 'mc','','medc','');
        [handlevecn{4},~,N{2},X{2},Y{2}] = dotdensity( kpos2, tP_oocv_h(id2_oocv), ...
            'dotEdgeColor', CLR, ...
            'dotFaceColor', CLR, ...
            'dotSize',MSoocv, ...
            'dotMarker', CLP, ...
            'meanLine', 'off', ...
            'medianLine', 'off');
        legvecn(4)=1;
        fN{2} = fid2_oocv(N{2});
    end

    if sum(~ind0)
        violin(P_oocv_h(~ind0),'x',5,'facecolor','k','edgecolor','none', 'facealpha', 0.3, 'bw',0.3);% 'mc','','medc','');
        [handlevecn{5},~,N{3},X{3},Y{3}] = dotdensity( 5 , P_oocv_h(~ind0), ...
            'dotEdgeColor', 'k', ...
            'dotFaceColor', 'k', ...
            'dotSize',MSoocv, ...
            'dotMarker', 'o', ...
            'meanLine', 'off', ...
            'medianLine', 'off');
        legvecn(5)=1;
        fN{3} = fid0_oocv(N{3});
    end
    N=cell2mat(fN');X=cell2mat(X');Y=cell2mat(Y');
else
    violin(P_oocv_h,'x',2,'facecolor','k','edgecolor','none', 'facealpha', 0.3, 'bw',0.3);% 'mc','','medc','');
    [handlevecn{5},~,N,X,Y] = dotdensity( 2 , P_oocv_h, ...
        'dotEdgeColor', 'k', ...
        'dotFaceColor', 'k', ...
        'dotSize',MSoocv, ...
        'dotMarker', 'o', ...
        'meanLine', 'off', ...
        'medianLine', 'off');
    legvecn(5)=1;
end

% Define textbox info data 
pss = cell(1,numel(P_oocv_h));psslen=0;
for i=1:numel(pss)
    expgroupi = 'not labeled';
    if labels_known
        if label_oocv_h(i) > 0
            expgroupi = handles.BinClass{h}.groupnames{1}; 
        elseif  label_oocv_h(i) < 0
            expgroupi = handles.BinClass{h}.groupnames{2}; 
        end
    end
    if sign(P_oocv_h(i)) > 0
        predgroupi = handles.BinClass{h}.groupnames{1}; 
    else
        predgroupi = handles.BinClass{h}.groupnames{2}; 
    end
    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g'], handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind}{i}, expgroupi, predgroupi, P_oocv_h(i));
    if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.cases       = handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind}(N);
figdata.x           = X;
figdata.y           = Y;
figdata.patterntext = pss(N);
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... 
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);

% Adjust y scaling depending on prob/rank estimates or decision scores
if GraphType >3 || handles.tglPercRank.Value || handles.tglP.Value
    probfx = zeroln;
    P_hx = P_h - probfx;
    handles.BinClass{h}.P_h = P_hx;
else
    probfx = 0;
    handles.BinClass{h}.P_h = P_h;
end

legvecn = logical(legvecn);
LegendStr = { sprintf('Discovery: %s', handles.BinClass{h}.groupnames{1}), ...
              sprintf('Discovery: %s', handles.BinClass{h}.groupnames{2}), ...
              sprintf('OOCV: %s', handles.BinClass{h}.groupnames{1}), ...
              sprintf('OOCV: %s', handles.BinClass{h}.groupnames{2}), ...
              'OOCV: unlabeled/not applicable'};
legendvec = LegendStr(legvecn);
handlevec = handlevecn(legvecn==1);
Hvec =[]; for i = 1:numel(handlevec), Hvec = [Hvec handlevec{i}]; end

% Plot binary class deviding line
xlims = numel(legendvec);
xLimits = get(handles.axes1,'XLim'); xLimitsVec = xLimits(1):xLimits(2); 
handles.axes1.XTick = 0.5:1:xlims+0.5; handles.axes1.YGrid = 'off'; handles.axes1.XGrid = 'on';
zeroline = ones(1,numel(xLimitsVec))*probfx;
plot(handles.axes1,xLimitsVec,zeroline,'k--','LineWidth',handles.ZeroLineWidth)
xlim(handles.axes1, [0.5 xlims+0.5]);
ylim(handles.axes1, 'auto');
handles.axes1.XTickLabel = [];

switch GraphType

    case {1,2,3}
        switch handles.params.TrainParam.SVM.prog
            case {'MikSVM','MKLRVM'}
                algostr = 'RVM probability';
            case 'LIBSVM'
                algostr = 'SVM score';
            case 'MVTRVR'
                algostr = 'RVR score';
            case 'MEXELM'
                algostr = 'ELM score';
            case 'LIBLIN'
                switch handles.params.TrainParam.SVM.LIBLIN.classifier
                    case {0,6}
                        algostr = 'LR';
                    otherwise
                        algostr = 'SVM';
                end
                switch handles.params.TrainParam.SVM.LIBLIN.b
                    case 0
                        algostr = [algostr ' Score'];
                    case 1
                        algostr = [algostr ' Probability'];
                end
            case 'matLRN'
                algostr = sprintf('matLearn [ %s ]',char(handles.params.TrainParam.SVM.matLRN.algo));
            otherwise
                algostr = [handles.params.TrainParam.SVM.prog 'score'];
        end
    case 4
        algostr = 'OOT-Probability';
    case 5
        algostr = 'Mean OOT-Probability (95%-CIs)';
    case 6
        algostr = 'Mean OOT-Probability (SD)';
end

if handles.tglPercRank.Value
    algostr = [algostr ' (%-ranks)'];
elseif handles.tglP.Value
    algostr = [algostr ' (prob)'];
end

hx(1) = xlabel('Samples'); 
set(hx(1), ...%'FontSize',handles.AxisLabelSize-2,
    'FontWeight',handles.AxisLabelWeight);
hx(2) = ylabel(['Binned ' algostr]); 
set(hx(2), ...%'FontSize',handles.AxisLabelSize-2, ...
    'FontWeight',handles.AxisLabelWeight);
handles.legend_classplot = legend(Hvec, legendvec, 'Location','Best','FontSize', 8,'LineWidth',1);%,'FontSize',handles.LegendFontSize); 
legend('boxon')
flg = 'off'; flg2='off';  
switch labels_known
    case 1
        flg='on';
        %% Extract contigency structure
        if subfl
            contigmat = ALLPARAM(label_oocv_h(ind0 & SubI), tP_oocv_h);
        else
            if isfield(handles.OOCV(oocvind).data,'BinResults')
                if iscell(handles.OOCV(oocvind).data.BinResults{1}.contingency)
                    contigmat = handles.OOCV(oocvind).data.BinResults{1}.contingency{h};
                else
                    contigmat = handles.OOCV(oocvind).data.BinResults{1}.contingency;
                end
            else
                contigmat = handles.OOCV(oocvind).data.MultiResults{1}.contingency{h};
            end
        end
        confmatrix = [[contigmat.TP contigmat.FN]; [contigmat.FP contigmat.TN]];

        %% Display Contigency data
        %if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
        handles = display_contigmat(handles, contigmat);

        %% Display contingency plot
        handles.h_contig = display_contigplot(handles, confmatrix, handles.BinClass{h}.groupnames);

        if sum(ind0)>2 %&& numel(unique(label_oocv_h(ind0)))>1 
            %% Display ROC
            if isfield(contigmat,'AUC')
                flg2 = 'on'; [handles.hroc, handles.hroc_random] = display_roc(handles, label_oocv_h (ind0 & SubI), tP_oocv_h);
            end
            %% Display pie charts
            [handles.h1pie, handles.h2pie] = display_piecharts(handles, [], contigmat, label_oocv_h(ind0 & SubI));   
        end
end
handles.pnRocCmds.Visible = flg2;
handles.pnPieCmds.Visible = flg2;
handles.pnContigCmds.Visible = flg;
handles.cmdMetricExport.Visible = flg;
handles.cmdExportAxes20.Visible = flg;
