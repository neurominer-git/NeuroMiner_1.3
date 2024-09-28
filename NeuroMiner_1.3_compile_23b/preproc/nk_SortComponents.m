function [TargParam_sorted, IN] = nk_SortComponents(matrix1, matrix2, IN)

% if isempty(IN) || ~isfield(IN,'actind')
%     for i=1:numel(TargParam)
%         if isfield(TargParam{i},'cmd') && strcmp(TargParam{i}.cmd,'reducedim')
%             IN.actind = i; break
%         end
%     end
% end
% if isempty(IN) || ~isfield(IN,'vecname')
%     IN.vecname = 'vec';
% end
% if isempty(IN) || ~isfield(IN,'mpp')
%     IN.mpp = 'mpp';
% end
%matrix1 = TemplParam{IN.actind}.(IN.mpp).(IN.vecname);
%matrix2 = TargParam{IN.actind}.(IN.mpp).(IN.vecname);

% Step 1: Calculate correlation between the columns of matrix1 and matrix2
correlations = corr(matrix1, matrix2);

% Step 2: Create a mapping between the columns of matrix1 and matrix2
numColumnsMatrix2 = size(matrix2, 2);

% Initialize mapping and flipping array with zeros
mapping = nan(numColumnsMatrix2, 1);
flipping = mapping;
%vec = 1:size(numColumnsMatrix1);

% For each column in matrix1, find the corresponding column index in matrix2
for colIndexMatrix2 = 1:numColumnsMatrix2
    
    if ~any(correlations(:))
        colIndexMatrix2=colIndexMatrix2-1; break, 
    end

    % Find the index of the column in matrix2 that has the highest correlation (positive or negative)
    [~, indexInMatrix1] = max(abs(correlations(:,colIndexMatrix2 )));

    % Store the corresponding column index from matrix2 in the mapping array
    mapping(colIndexMatrix2) = indexInMatrix1;
    flipping(colIndexMatrix2) = sign(correlations(indexInMatrix1, colIndexMatrix2));

    % If the chosen column in matrix2 is already assigned, find the next best unassigned column
    correlations(indexInMatrix1, :) = nan; % Set correlation to Nan to exclude it from consideration
end
idxnan = isnan(mapping);
mapping(idxnan) = [];
flipping(idxnan) =[];
% Step 3: Reorder the columns of matrix2 based on the mapping
IN.templmatrixsize = size(matrix1);
IN.correlations = correlations;
IN.mapping = mapping;
IN.flipping = flipping;
try
    TargParam_sorted = nan( IN.templmatrixsize );
    TargParam_sorted(:,mapping) = matrix2(:,1:colIndexMatrix2);
catch
    fprintf('problem')
end