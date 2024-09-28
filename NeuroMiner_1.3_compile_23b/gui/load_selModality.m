function load_selModality(handles)

if strcmp(handles.popupmenu1.String{handles.popupmenu1.Value},'Multi-group classifier')
    multiflag = true;
else
    multiflag = false;
end

popuplist=[];
if ~multiflag
    for i=1:size(handles.visdata,1)
        popuplist{i} = sprintf('Modality %g: %s', i, handles.NM.datadescriptor{handles.visdata{i,handles.curlabel}.params.varind}.desc);
    end
    handles.selModality.String = popuplist;
    handles.selModality.Value = handles.curmodal;
    if handles.curmodal<=size(handles.visdata,1)
        popuplist=[];
        v = handles.visdata{handles.curmodal,handles.curlabel};
        if isfield(v,'CVRnorm')
            CVRnorm_opts = {'(SD-based) ', '(SEM-based) '};
            CVRnormStr = CVRnorm_opts{v.CVRnorm};
        else
            CVRnormStr = '';
        end
        if isfield(v,'MEAN'),                      popuplist{1} = 'Feature weights [Overall Mean (StErr)]';                          end
        if isfield(v,'MEAN_CV2'),                  popuplist{end+1} = 'Feature weights [CV2 Mean (StErr)]';                          end
        if isfield(v,'CVRatio'),                   popuplist{end+1} = sprintf('CVR %sof feature weights [Overall Mean]',CVRnormStr); end
        if isfield(v,'CVRatio_CV2'),               popuplist{end+1} = sprintf('CVR %sof feature weights [CV2 Mean]', CVRnormStr);    end
        if isfield(v,'FeatProb'),                  popuplist{end+1} = 'Feature selection probability [Overall Mean]';                end
        if isfield(v,'Prob_CV2'),                  popuplist{end+1} = 'Probability of feature reliability (95%-CI) [CV2 Mean]';      end
	    if isfield(v,'SignBased_CV2'),             popuplist{end+1} = 'Sign-based consistency';  								     end
	    if isfield(v,'SignBased_CV2_z'),           popuplist{end+1} = 'Sign-based consistency (Z score)';  						     end
	    if isfield(v,'SignBased_CV2_p_uncorr'),    popuplist{end+1} = 'Sign-based consistency -log10(P value)';  				     end
	    if isfield(v,'SignBased_CV2_p_fdr'),       popuplist{end+1} = 'Sign-based consistency -log10(P value, FDR)';  			     end
        if isfield(v,'Spearman_CV2'),              popuplist{end+1} = 'Spearman correlation [CV2 Mean]';                             end
        if isfield(v,'Pearson_CV2'),               popuplist{end+1} = 'Pearson correlation [CV2 Mean]';                              end
        if isfield(v,'Spearman_CV2_p_uncorr'),     popuplist{end+1} = 'Spearman correlation -log10(P value) [CV2 Mean]';             end
        if isfield(v,'Pearson_CV2_p_uncorr'),      popuplist{end+1} = 'Pearson correlation -log10(P value) [CV2 Mean]';              end
        if isfield(v,'Spearman_CV2_p_fdr'),        popuplist{end+1} = 'Spearman correlation -log10(P value, FDR) [CV2 Mean]';        end
        if isfield(v,'Pearson_CV2_p_fdr'),         popuplist{end+1} = 'Pearson correlation -log10(P value, FDR) [CV2 Mean]';         end
        if isfield(v,'PermProb_CV2'),              popuplist{end+1} = 'Permutation-based -log10(P value) [CV2 Mean]';                end
        if isfield(v,'PermProb_CV2_FDR_PVAL'),     popuplist{end+1} = 'Permutation-based -log10(P value, FDR) [CV2 Mean]';           end
        if isfield(v,'PermZ_CV2'),                 popuplist{end+1} = 'Permutation-based Z Score [CV2 Mean]';                        end
	    if isfield(v,'Analytical_p'),     		   popuplist{end+1} = 'Analytical -log10(P Value) for Linear SVM [CV2 Mean]';        end
	    if isfield(v,'Analyitcal_p_fdr'),     	   popuplist{end+1} = 'Analytical -log10(P Value, FDR) for Linear SVM [CV2 Mean]';   end
        if isfield(v,'PermModel_Eval_Global'),     popuplist{end+1} = 'Model P value histogram';                                     end
        if isfield(v,'ExtraL')
            for i=1:numel(v.ExtraL)
                popuplist{end+1} = sprintf('Generalization analysis for extra label ''%s'' [#%g]',v.ExtraL(i).LABEL_NAME, i);
            end
        end
        if isfield(v,'CorrMat_CV2'),               popuplist{end+1} = 'Correlation matrix';                                        end
        if isfield(v,'CorrMat_CV2_p_uncorr'),      popuplist{end+1} = 'Correlation matrix (P value)';                              end
        if isfield(v,'CorrMat_CV2_p_fdr'),         popuplist{end+1} = 'Correlation matrix (P value, FDR)';                       end
        if isfield(v,'CorrMat_CV2'),               popuplist{end+1} = 'Network plot correlation matrix';                           end
        if isfield(v,'CorrMat_CV2_p_uncorr'),      popuplist{end+1} = 'Network plot correlation matrix (P value)';                 end
        if isfield(v,'CorrMat_CV2_p_fdr'),         popuplist{end+1} = 'Network plot correlation matrix (P value, FDR)';            end
        handles.selVisMeas.String = popuplist; 
        VisOnFl = 'on';
        VisElFl = 'on';
    else
        VisOnFl = 'off';
        VisElFl = 'off';
        handles.curmodal = size(handles.visdata,1);
    end
    
else
    v = handles.visdata{handles.curmodal,handles.curlabel};
    if isfield(v,'ObsModel_Eval_Global_Multi')
        handles.selVisMeas.Enable = 'on'; 
        popuplist{1} = 'Model P value histogram [Multi-group]';
        for i=1:numel(v.ObsModel_Eval_Global_Multi_Bin)
            popuplist{end+1} = sprintf('Model P value histogram [Class %s vs. Rest]',handles.NM.groupnames{i});    
        end
        handles.selVisMeas.String = popuplist;
        VisOnFl = 'on';
        VisElFl = 'off';
    else
        VisOnFl = 'on';
        VisElFl = 'on';
    end
end
% Toggle enabled status of visualization controls
handles.selVisMeas.Enable = VisOnFl; 
switch handles.params.TrainParam.FUSION.flag
    case 3
        handles.selModality.Enable = 'off';
    otherwise
        handles.selModality.Enable = VisElFl;
end
handles.selPager.Enable = VisElFl;
handles.tglSortFeat.Enable = VisElFl;
handles.cmdExportFeats.Enable = VisElFl;