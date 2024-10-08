function [f, ax, ax2] = barShapley(casenum, MLIcont, feats, refdata, sortfl, uiaxes, uiaxes2, col1, col2) %newdata

if exist('uiaxes','var')
    ax = uiaxes;
    f=0;
else
    f=figure;
    f.Position = [100 100 750 750];
    ax = gca;
end

hold(ax,'on');
nF = numel(feats);
y = MLIcont.ShapleyValues(casenum,:);

if sortfl
    [~,idx] = sort(y);
else
    idx = 1:size(y,2);
end


b1 = bar(ax, y(idx),'FaceColor',rgb('DodgerBlue'),'EdgeColor',rgb("SlateGrey"));

legend(ax, b1, 'Shapley Values');



%ax=gca;
ax.XTick=1:nF;

feats = regexprep(feats,'_', '\\_');
ax.XTickLabel=feats(idx);
ax.YAxis.Label.String = 'Shapley value';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 10;
ax.XAxis.Label.String = 'Features';
ax.Box="on";
ax.XTickLabelRotation=15;
if ~isnan(max(y))
    ax.YLim=[min(min(y),0)*1.2,max(y)*1.2];
end
if exist('refdata','var') && ~isempty(refdata)
    if exist('uiaxes2','var')
        ax2 = uiaxes2;
        % resize to size of other plot
        ax2.InnerPosition(1) = ax.InnerPosition(1);
        ax2.InnerPosition(3) = ax.InnerPosition(3);
    else
        ax.Position(4) = .50;
        pos = ax.Position;
        pos([2 4]) = [0.75 0.2];
        ax2 = axes('Position', pos);
    end
    if exist('origdata','var') && ~isempty(origdata) % for CARE project
        refdata = [origdata(casenum,:), refdata];
        centiles = nk_ComputePercentiles(refdata,origdata(casenum,:),'inverse');
    else
        centiles = nk_ComputePercentiles(refdata,refdata(casenum,:),'inverse');
        grouplabel = evalin('base', 'NM.label');
        group1_idx = grouplabel == 1;
        group2_idx = grouplabel == 2; 
        group1 = refdata(group1_idx,:);
        
        for i = 1:size(group1,1)
            centiles_i = nk_ComputePercentiles(refdata, group1(i,:), 'inverse');
            if i == 1
                centiles_group1 = centiles_i;
            else 
                centiles_group1 = [centiles_group1; centiles_i];
            end
        end
        av_centiles_group1 = mean(centiles_group1);

        group2 = refdata(group2_idx,:);
        for i = 1:size(group2,1)
            centiles_i = nk_ComputePercentiles(refdata, group2(i,:), 'inverse');
            if i == 1
                centiles_group2 = centiles_i;
            else 
                centiles_group2 = [centiles_group2; centiles_i];
            end
        end
        av_centiles_group2 = mean(centiles_group2);
        
    end
    
    bar(ax2, centiles(idx),'FaceColor',rgb('SlateGray'),'EdgeColor',rgb("Black"));
    hold(ax2,'on');
    plot(ax2, av_centiles_group1(idx), 'LineWidth', 3, 'color', col1);%,'FaceColor',rgb('SlateGray'),'EdgeColor',rgb("Black"));
    plot(ax2, av_centiles_group2(idx), 'LineWidth', 3, 'color', col2);
    hold(ax2, 'off');

    groupnames = evalin('base', 'NM.groupnames');
    legend(ax2, "subject's percentiles" , sprintf('mean percentiles of group %s', char(groupnames(1))), ...
        sprintf('mean percentiles of group %s', char(groupnames(2))));

    ax2.YAxis.Label.String = 'Percentile rank [%]';
    ax2.YAxis.Label.FontWeight = 'bold';
    ax2.YAxis.Label.FontSize = 10;
    ax2.XTick=1:numel(feats);
    ax2.XLim=[-0.2 numel(feats)+1.2];
    ax2.Box ="on";
    ax2.XTickLabel=[];
    if ~exist('uiaxes2','var')
        g = grouplabel(casenum);
        case_groupname = char(groupnames(g));
        ax2.Title.String = sprintf('Shapley values of subject #%g, Group: %s', casenum, case_groupname);
        ax2.Title.FontWeight = 'bold';
        ax2.Title.FontSize = 14;
    end
    ax2.XTickLabel=feats(idx);
    ax2.XTickLabelRotation=15;
    ax2.YLim=[ 0 100 ];
else
    if ~exist('uiaxes','var')
        ax.Title.String = sprintf('Shapley values of subject #%g', casenum);
    end
    ax.Title.FontWeight = 'bold';
    ax.Title.FontSize = 14;
end
