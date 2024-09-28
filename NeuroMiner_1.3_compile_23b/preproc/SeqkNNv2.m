function imputed_data = SeqkNNv2(data, K, traindata)
% SeqKNN: Sequential KNN imputation method
% This function estimates missing values sequentially from the gene that has
% least missing rate in microarray data, using weighted mean of k nearest neighbors.
%
% <Usage>
% imputed_data=SeqKNN(data, k, traindata);
%
% <Arguments>
% data: matrix or dataframe, 1 row corresponds to 1 gene, 1 column to 1
% sample, colnames and rownames can be used
% K: number of nearest neighbors
% traindata: (optional) training data to use when all rows in data have missing values
%
% <Details>
% SeqKNN separates the dataset into incomplete and complete sets that have
% or have not missing values, respectively. 
% The genes in incomplete set are imputed by the order of missing rate. Missing value
% is filled by the weighted mean value of corresponding column of the nearest neighbor genes in
% complete set. Once all missing values in a gene are imputed, the imputed gene is moved into the
% complete set and used for the imputation of the rest of genes in incomplete set.

imputed_data = zeros(size(data));
complete = [];
incomplete = [];
missing = [];
com_ind = [];
incom_ind = [];
[rows, cols] = size(data);

for i = 1:rows
    if ~any(isnan(data(i,:)))
        complete = [complete; data(i,:)];
        com_ind = [com_ind, i];
    else
        incomplete = [incomplete; data(i,:)];
        incom_ind = [incom_ind, i];
        missing = [missing, sum(isnan(data(i,:)))];
    end
end

% Fallback imputation for completely missing rows if no complete rows are available
if isempty(complete) || exist('traindata', 'var')
    if exist('traindata', 'var') && ~sum(isnan(traindata(:)))
        complete = traindata;
    else
        % Check if there is at least one complete row
        % Use column-wise mean of the current dataset
        meanValues = repmat(nm_nanmean(data), rows, 1);
        % Replace NaNs with column-wise means
        idx = isnan(data);
        complete = data;
        complete(idx) = meanValues(idx);
    end
end
imputed_data(com_ind,:) = data(com_ind,:);
[~, missing_order] = sort(missing);
completedata = complete;

% Part 2: Impute missing values
for j = 1:size(incomplete,1)
    fprintf('.');
    dist = [];
    cgen = size(completedata, 1);
    for i = 1:cgen
        % Calculate distance, handling NaNs appropriately
        dist(i) = nansum((incomplete(missing_order(j),:) - completedata(i,:)).^2);
    end
    [dist, pos] = sort(dist);
    pos = pos(1:K);
    dist = dist(1:K);
    dist(dist == 0) = realmin; % Avoid division by zero
    weights = (1./dist) ./ sum(1./dist);
    for g = 1:cols
        if isnan(incomplete(missing_order(j),g))
            incomplete(missing_order(j),g) = sum(weights .* completedata(pos,g)');
        end
    end
    completedata = [completedata; incomplete(missing_order(j),:)];
end

imputed_data(incom_ind,:) = incomplete;
end