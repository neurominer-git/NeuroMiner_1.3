function NM = update_NMstruct(NM)
% function that updates an older NM structure to be compatible with the new
% release

% add label-field to analysis if necessary

% create default label struct
label.label = NM.label; 
label.modeflag = NM.modeflag; 
label.altlabelflag = 0; 

for i=1:length(NM.analysis)
    if ~isfield(NM.analysis{i}.params, 'label')
        NM.analysis{1,i}.params.label = label;
    end
    if isfield(NM.analysis{i},'visdata') && iscell(NM.analysis{i}.visdata{1})
        for j=1:numel(NM.analysis{i}.visdata)
            for k=1:size(NM.analysis{i}.visdata{j})
                visdata{j,k} = NM.analysis{i}.visdata{j}{k};
            end
        end
        NM.analysis{i}.visdata = visdata;
    end
end

