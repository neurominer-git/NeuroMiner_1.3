function [sY, sX, sYX, IN] = nk_PLS(Y, X, IN)
% Computes PLS model from experimental predictors X and predictands Y
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 10/2023
if ~isfield(IN,'trained'),IN.trained=0;end
if ~isfield(IN,'algostr'),IN.algostr='pls';end

% Get predictor matrix data (X)
if exist('X','var') && ~isempty(X)
    Xd = X;
else
    Xd = IN.DR.PLS.V;
end

if ~IN.trained % No LV projection parameters found => compute mpp
    
    % Step 1
    % Mean centering / standardization of predictands (Y) and predictor
    % matrix (Xd) depending on the deviation mapping methodology
    switch IN.algostr
        case 'pls'
            IN.mpp.mY.method = 'mean-centering';
            IN.mpp.mX.method = 'mean-centering';
        case 'spls'
            IN.mpp.mY.method = 'standardization using mean';
            IN.mpp.mX.method = 'standardization using mean';
    end
    [mY, IN.mpp.mY] = nk_PerfStandardizeObj(Y,IN.mpp.mY);
    [mXd, IN.mpp.mX] = nk_PerfStandardizeObj(Xd,IN.mpp.mX);
    
    % Step 2 (optional step)
    % If mXd contains NaNs, we have to implement an imputation step here:
    % Imputation
    if sum(isnan(mXd(:)))
        mXd = SeqkNN(mXd, 5);
        % if there are still missings use median replacement
        if sum(isnan(mXd(:)))
             IN.mn = median(mXd,'omitnan');
             for i = 1:size(mXd,2)
                 idxnan = isnan(mXd(:,i));
                 mXd(idxnan,i) = IN.mn(i);
             end
        end
        IN.tX = mXd;
    end
    
    % Step 3 (learn latent predictor space: Y = U_1 * d_1 * V_1 + U_2 * d_2 * V_2 + ... U_n * d_n * V_n)
    switch IN.algostr
        case 'pls'
            % Build covariance matrix S
            S=mXd'*mY;
            % Decompose S into LVpairs
            [IN.mpp.u, IN.mpp.s, IN.mpp.v] = svd(S',0);
            IN.mpp.d = diag(IN.mpp.s)';
        case 'spls'
            % Perform sparse PLS
            nD = size(mXd,2);
            IN.mpp.u = zeros(size(mY,2),nD);
            IN.mpp.v = zeros(nD,nD);
            IN.mpp.C = [];
            for i=1:nD
                % Compute SPLS level i
                [IN.mpp.u(:,i), IN.mpp.v(:,i), IN.mpp.C] = spls(mY, mXd, IN.cu, IN.cv, IN.mpp.C);
                % Covariance matrix deflation
                IN.mpp.C = IN.mpp.C - (IN.mpp.C * IN.mpp.v(:,i)) * IN.mpp.v(:,i)';
            end
    end
    IN.trained = true;
else
    % Mean-center the Y (predictand) matrix using the parameters learned in the training
    % data
    mY = nk_PerfStandardizeObj(Y,IN.mpp.mY);
    % Compute U scores by applying U transform to the mean-centered test
    % data
    if ~isempty(Xd)
        % Mean-center the X (predictor) matrix using the parameters learned in the training
        % data
        mXd = nk_PerfStandardizeObj(Xd,IN.mpp.mX);
        % If mXd contains NaNs, we have to implement an imputation step here:
        % Imputation
        if sum(isnan(mXd(:)))
            % Perform sequential kNN imputation
            mXd = SeqkNN(mXd, 5, IN.tX);
            % if there are still missings use median replacement
            if sum(isnan(mXd(:)))
                 if ~isfield(IN,'mn')
                     mn = median(IN.tX,'omitnan');
                 else
                     mn=IN.mn;
                 end
                 for i = 1:size(mXd,2)
                     idxnan = isnan(mXd(:,i));
                     mXd(idxnan,i) = mn(i);
                 end
            end
        end
    end
end

% Final step

% Compute U scores by applying U transform to the mean-centered
% predictand matrix (Y)
sY = mY*IN.mpp.u;
% Now compute V scores by applying the learned V transform to
% the mean-centered predictor matrix (X).
sX = mXd*IN.mpp.v;
% Compute UV scores for the predictor data
sYX = mXd * (IN.mpp.u * IN.mpp.v')' + IN.mpp.mY.meanY;
