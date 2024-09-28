function R = nk_SVMFeatRank(Y, label, SVM, ModelOnly)
global MODEFL

if ~isfield(SVM,'modeflag')
    modeflag = MODEFL;
else
    modeflag = SVM.modeflag;
end

CMDSTR = nk_DefineCmdStr(SVM, modeflag);
[~,~,EVALFUNC,~] = nk_GetScaleYAxisLabel(SVM.evalfunc);
if strcmp(SVM.prog,'LIBSVM') 
    [LIBSVMTRAIN, LIBSVMPREDICT] = nk_DefineLIBSVMfun(SVM);
end
if ~exist('ModelOnly','var') || isempty(ModelOnly)
    ModelOnly = 1;
end

if strcmp(modeflag,'classification')
    uLabels = unique(label);
    nu = numel(uLabels);
else
    nu = 1;
end

P = nk_CreateSVMParamArray(SVM);
Ps = allcomb(P.Params,'matlab');
[nP, ~] = size(Ps);

R = zeros(size(Y,2),nP);
if size(label,1) ~= size(Y,1), label = label'; end
RStr = ''; SlackStr = ''; 

if strcmp(SVM.prog,'LIBLIN'), Y = sparse(Y); end
    
for i = 1:nP
    
    %% Train SVM model
    switch SVM.prog
        case 'LIBSVM'
            %% Give verbose info regarding training parameters
            % Slack parameter (either C or Nu)
            if strcmp(P.Params_desc{1},'SlackParam')
                Pt.SlackParam = num2str(Ps(i,1),'%1.5f');
                SlackStr = sprintf('Slack parameter: %s', Pt.SlackParam);
            elseif strcmp(P.Params_desc{1},'NuCParam')
                Pt.NuCParam = num2str(Ps(i,1),'%1.5f');
                SlackStr = sprintf('NuC parameter: %s', Pt.NuCParam);
            end
            % Regression accuracy parameter (either Epsilon or Nu)
            if strcmp(SVM.modeflag,'regression')
                if strcmp(P.Params_desc{2},'NuParam')
                    Pt.NuParam = num2str(Ps(i,2),'%1.5f'); 
                    RStr = sprintf(', Nu parameter (SVR): %s', Pt.NuParam);
                elseif strcmp(P.Params_desc{2},'EpsParam')
                    Pt.EpsParam = num2str(Ps(i,2),'%1.5f'); 
                    RStr = sprintf(', Epsilon parameter (SVR): %s', Pt.EpsParam);
                end
            end
            fprintf('\nTraining linear LIBSVM model with %s%s%s', SlackStr, RStr);
            [~, model] = nk_GetParamSVM(Y, label, CMDSTR, SVM, LIBSVMTRAIN, LIBSVMPREDICT, EVALFUNC, Pt, ModelOnly);
            fprintf(' => Model complexity = %g%%',model.totalSV*100/numel(label));
    
            %% Print performance
            %fprintf(' ==> %s = %g', EVALFUNC, param.val);

            %% Get weight vector
            if nu < 3 % binary
                % Compute primal variable w
                R(:,i) = model.SVs' * model.sv_coef;

            else % Multi-class with one-vs-one decomposition (NOT WORKING!)
                k = nu*(nu-1)/2;
                Ri = zeros(size(Y,2), k); cnt = 1;
                for j = 1:k
                    %+-+-+--------------------+
                    %|1|1|                    |
                    %|v|v|  SVs from class 1  |
                    %|2|3|                    |
                    %+-+-+--------------------+
                    %|1|2|                    |
                    %|v|v|  SVs from class 2  |
                    %|2|3|                    |
                    %+-+-+--------------------+
                    %|1|2|                    |
                    %|v|v|  SVs from class 3  |
                    %|3|3|                    |
                    %+-+-+--------------------+

                    % Get SV index of positive class
                    if j == 1
                        indP = 1:model.nSV(j);
                    else
                        indP = sum(model.nSV(1:j-1))+1 : sum(model.nSV(1:j)) ;
                    end

                    for jl = j+1:k

                        % Get SV index of negative class
                        indN = sum(model.nSV(1:jl-1))+1 : sum(model.nSV(1:jl));
                        coef = [model.sv_coef(indP,2); model.sv_coef(indN,1)];
                        SVs = [model.SVs(indP,:); model.SVs(indN,:)];
                        Ri(:,cnt) = SVs'*coef;
                        cnt = cnt+1;

                    end

                end
                % Now compute average weight vector
                R(:,i) = mean(Ri,2);

            end
            
        case 'LIBLIN'
            EpsParam = ''; Tolerance = ' -e 0.01'; SlackParam = ' -c 1';
            idx = find(strcmp(P.Params_desc,'SlackParam')); 
            if idx, SlackParam = sprintf( ' -c %1.4f', Ps(i,idx)); end
            idx = find(strcmp(P.Params_desc,'EpsParam')); 
            if idx, EpsParam = sprintf( ' -p %1.4f', Ps(i,idx)); end
            idx = find(strcmp(P.Params_desc,'Tol')); 
            if idx, Tolerance = sprintf( ' -e %1.4f', Ps(i,idx)); end
            cmdstr = [  SlackParam EpsParam Tolerance CMDSTR.notolmodel CMDSTR.quiet];
            if nu>2
                pattern = '-s\s[1-9]';
                cmdstr = regexprep(cmdstr,pattern,'-s 4');
                fprintf('\Feature weighting using multi-class LIBLINEAR model using Crammer-Singer algorithm: %s', cmdstr);
            else
                fprintf('\nFeature weighting using LIBLINEAR model with command string:%s', cmdstr);
            end
            model = train_liblin244(label,Y,cmdstr);
            if nu < 3 % binary
                if contains(cmdstr,'-B 1')
                    R(:,i) = model.w(1:end-1); 
                else
                    R(:,i) = model.w;
                end
            else
                if contains(cmdstr,'-B 1')
                    R(:,i) = sum(model.w(:,1:end-1)); 
                else
                    R(:,i) = sum(model.w); % multiclass
                end
            end
            nonzerofeat = sum(R(:,i) ~= 0); 
            fprintf(' ==> %g (%1.2f%%) Non-zero features.', nonzerofeat, nonzerofeat*100/model.nr_feature);
    end
    
end

end