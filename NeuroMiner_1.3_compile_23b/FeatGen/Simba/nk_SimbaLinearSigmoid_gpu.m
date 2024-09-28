function [alphas, margin_list] = nk_SimbaLinearSigmoid_gpu(X_train,  Y_train, extra_param)

% alphas = simbaOneTime(X_train,  Y_train, extra_param);
%
% This function runs simba with linear utility function one time (i.e. with one starting point).
%
% input: X_train(i,j) is the value of feature j in training instance i
%        Y_train(i) is the label of training instance i
%        extra_param is a struct that may contain the following parameters for the
%        algorithm: 
%                  max_iter: the number of pathes on all training data. (default is 1)
%                  blocksize: the algorithm recaculate the distances after updating the features weights
%                             usig "blocksize" samples. (default is 1)
%                  verbose: 1 for verbose, 0 otherwise (default is 0)
%
% output: alphas(j) is the weight of the j's feature. higher is better.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written by Amir Navot & Ran Gilad-Bachrach                   %%
%% Date: April 1, 2004                                           %%
%% Last update: August 2nd, 2004 by Ran Gilad-Bachrach           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modified version run on CUDA enabled GPUs to save time
% ======================================================
% Nikolaos Koutsouleris 11/2009

if isfield(extra_param, 'verbose'), verbose = extra_param.verbose; else, verbose = 0;  end
if isfield(extra_param, 'blocksize'), blocksize = extra_param.blocksize; else, blocksize = 1; 
    if verbose,disp(['blocksize = ' num2str(blocksize)]); end; end
if isfield(extra_param, 'max_iter'), max_iter = extra_param.max_iter; else, max_iter = 1; 
    if verbose,disp(['max_iter = ' num2str(max_iter)]); end; end
if isfield(extra_param, 'iter_abort_crit'), iter_abort_crit = extra_param.iter_abort_crit; else, iter_abort_crit = 100; 
    if verbose,disp(['iter_abort_crit = ' num2str(iter_abort_crit)]); end; end
if isfield(extra_param,'beta'), sigmoidfl = true ;else sigmoidfl=false; end

feat_num = size(X_train,2);
N = size(X_train,1);
if blocksize > N,
    error('blocksize > N');
end

%alphas = ones(1,feat_num);
featind = ones(1,feat_num);
g_alphas = ones(1,feat_num, GPUsingle);
% Move / create variables on the GPU
g_X_train = GPUsingle(X_train);
g_Y_train = GPUsingle(Y_train);
g_X_train_square = g_X_train.^2;
labelsSet = unique(Y_train);
margin_list = zeros(max_iter,1);

for iter = 1:max_iter,
    tic;
    perm = randperm(N);

    for bi = 0:blocksize:(N-blocksize),

        if(verbose && ~mod(bi,10))
            disp(bi);
        end;

        xi = (bi+1):(bi+blocksize);
        xs = perm(xi);
        
        [pplus , pminus, delta_plus, delta_minus, margins] = ...
            nk_SimbaHelperGpu(g_X_train, g_X_train_square, g_Y_train, Y_train, g_alphas, xs, labelsSet);
   
        J = find(delta_plus~=0 & delta_minus~=0);
        if (length(J)~=length(delta_minus))
            warning('distance to nearhit or nearmiss is 0');
            if isempty(J),
                continue;
            end
        end
        xs = xs(J);
        pplus = pplus(J);
        pminus = pminus(J);
        delta_plus = delta_plus(J);
        delta_minus = delta_minus(J);
        delta_plus = delta_plus(:);
        delta_minus = delta_minus(:);
        vecs1 = (X_train(xs,:) - X_train(pplus,:)).^2;
        vecs2 = (X_train(xs,:) - X_train(pminus,:)).^2;
        vecs1 = vecs1 ./ delta_plus(:, featind);
        vecs2 = vecs2 ./ delta_minus(:, featind);
        
        if ~sigmoidfl
            dw = (1/blocksize)*(sum(vecs2,1) - sum(vecs1,1));
        else
            ebx = exp(-beta*(delta_minus-delta_plus)/2);
            dl = beta * ebx ./ ((1+ebx).^2);
            Dl = dl(:, featind);
            Dlw = (vecs2 - vecs1) .* Dl;
            dw = (1/blocksize)*(sum(Dlw,1));
        end
        
        g_alphas = g_alphas .* GPUsingle(1 + dw);

    end
    if bi < (N-blocksize),
        
        xi = (bi+blocksize+1):N;
        xs = perm(xi);
        
        [pplus , pminus, delta_plus, delta_minus, margins] = ...
            nk_SimbaHelperGpu(g_X_train, g_X_train_square, g_Y_train, Y_train, g_alphas, xs, labelsSet);
    
        delta_plus = delta_plus(:);
        delta_minus = delta_minus(:);
        vecs1 = (X_train(xs,:) - X_train(pplus,:)).^2;
        vecs2 = (X_train(xs,:) - X_train(pminus,:)).^2;
        vecs1 = vecs1 ./ delta_plus(:, featind);
        vecs2 = vecs2 ./ delta_minus(:, featind);
        if ~sigmoidfl
            dw = (1/blocksize)*(sum(vecs2,1) - sum(vecs1,1));
        else
            ebx = exp(-beta*(delta_minus-delta_plus)/2);
            dl = beta * ebx ./ ((1+ebx).^2);
            Dl = dl(:, featind);
            Dlw = (vecs2 - vecs1) .* Dl;
            dw = (1/blocksize)*(sum(Dlw,1));
        end
        
        g_alphas = g_alphas .* GPUsingle(1 + dw);

    end;
    tElapsed = toc;
    margin_list(iter) = min(margins);
    if perc_poschange == iter_abort_crit, break; end
    fprintf('\nIter %g/%g:\tfeature stability index = %1.3f%% (# unstable features = %g).\tDone in %1.2f sec.', ...
        iter , max_iter, perc_poschange, abs_poschange, tElapsed)    
end

alphas = single(g_alphas);

