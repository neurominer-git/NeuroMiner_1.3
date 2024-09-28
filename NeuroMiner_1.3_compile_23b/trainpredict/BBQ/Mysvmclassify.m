function outclass = Mysvmclassify(svmStruct,sample, varargin)
%SVMCLASSIFY classifies data using a support vector machine
%
%   GROUP = SVMCLASSIFY(SVMStruct,SAMPLE) classifies each row of the data
%   in SAMPLE using the information in a support vector machine classifier
%   structure SVMStruct created using SVMTRAIN. SAMPLE must have the same
%   number of columns as the data used to train the classifier in SVMTRAIN.
%   GROUP indicates the group to which each row of SAMPLE is assigned.
%
%   SVMCLASSIFY(...,'SHOWPLOT',true) plots the sample data on the figure
%   created using the SHOWPLOT option in SVMTRAIN.
%
%   Example:
%       % Load the data and select features for classification
%       load fisheriris
%       data = [meas(:,1), meas(:,2)];
%       % Extract the Setosa class
%       groups = ismember(species,'setosa');
%       % Randomly select training and test sets
%       [train, test] = crossvalind('holdOut',groups);
%       cp = classperf(groups);
%       % Use a linear support vector machine classifier
%       svmStruct = svmtrain(data(train,:),groups(train),'showplot',true);
%       classes = svmclassify(svmStruct,data(test,:),'showplot',true);
%       % See how well the classifier performed
%       classperf(cp,classes,test);
%       cp.CorrectRate
%
%   See also CLASSIFY, CLASSPERF, CROSSVALIND, KNNCLASSIFY, QUADPROG,
%   SVMSMOSET, SVMTRAIN.

%   Copyright 2004-2008 The MathWorks, Inc.
%   $Revision: 1.1.12.11 $  $Date: 2008/06/16 16:32:46 $

%   References:
%
%     [1] Cristianini, N., Shawe-Taylor, J An Introduction to Support
%         Vector Machines, Cambridge University Press, Cambridge, UK. 2000.
%         http://www.support-vector.net
%     [2] Kecman, V, Learning and Soft Computing,
%         MIT Press, Cambridge, MA. 2001.
%     [3] Suykens, J.A.K., Van Gestel, T., De Brabanter, J., De Moor, B.,
%         Vandewalle, J., Least Squares Support Vector Machines,
%         World Scientific, Singapore, 2002.


% set defaults


plotflag = false;

% check inputs
bioinfochecknargin(nargin,2,mfilename);

% deal with struct input case
if ~isstruct(svmStruct)
    error('Bioinfo:svmclassify:TwoInputsNoStruct',...
        'The first input should be a struct generated by SVMTRAIN.');
end


% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error('Bioinfo:svmclassify:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'showplot','-compilerhelper'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:svmclassify:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:svmclassify:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % plotflag
                    plotflag = opttf(pval,okargs{k},mfilename);
                case 2 % help the compiler find required function handles by including svmtrain
                    svmtrain(eye(2),[1 0]);
            end
        end
    end
end

groupnames = svmStruct.GroupNames;

% check group is a vector -- though char input is special...
if ~isvector(groupnames) && ~ischar(groupnames)
    error('Bioinfo:svmclassify:GroupNotVector',...
        'Group must be a vector.');
end

% grp2idx sorts a numeric grouping var ascending, and a string grouping
% var by order of first occurrence

[g,groupString] = grp2idx(groupnames);  

% do the classification
if ~isempty(sample)
    % shift and scale the data if necessary:
    sampleOrig = sample;
    if ~isempty(svmStruct.ScaleData)
        for c = 1:size(sample, 2)
            sample(:,c) = svmStruct.ScaleData.scaleFactor(c) * ...
                (sample(:,c) +  svmStruct.ScaleData.shift(c));
        end
    end

    try
%         Edited By Mahdi
        classified = Mysvmdecision(sample,svmStruct);
%         End Of Edition By Mahdi
    catch theException
        error('Bioinfo:svmclassify:ClassifyFailed',...
            'An error was encountered during classification.\n%s',theException.message);
    end
    if plotflag

        if isempty(svmStruct.FigureHandles)
            warning('Bioinfo:svmclassify:NoTrainingFigure',...
                'No figure was created during training so no ShowPlot option does not work.');

        else
            try
                hAxis = svmStruct.FigureHandles{1};
                hLines = svmStruct.FigureHandles{2};
                hSV = svmStruct.FigureHandles{3};
                % unscale the data for plotting purposes
                [hAxis,hClassLines] = svmplotdata(sampleOrig,classified,hAxis); 
                trainingString = strcat(cellstr(groupString),' (training)');
                sampleString = strcat(cellstr(groupString),' (classified)');
                legend([hLines(1),hClassLines(1),hLines(2),hClassLines(2),hSV],...
                    {trainingString{1},sampleString{1},...
                    trainingString{2},sampleString{2},'Support Vectors'});
            catch theException
                warning('Bioinfo:svmclassify:DisplayFailed',...
                    'An error was encountered during plotting.\n%s',theException.message);
            end
        end
    end
    
%     Commented by Mahdi
    %     classified(classified == -1) = 2;
%     End of Edition by Mahdi

    outclass = classified;
    unClassified = isnan(outclass);
    % if there ara unclassified points (NaNs) deal with the situation...
    if any(unClassified)
        warning('Bioinfo:svmclassify:UnclassifiedData',...
            'Some samples could not be classified. This is probably caused by NaN values in the data.')
        numGroups = numel(groupString);
        % for numeric groups we use NaN for unclassifiable, for string
        % group names use empty string.
        if isnumeric(groupnames) || islogical(groupnames)
            groupString{end+1} = 'NaN';
        else
            groupString{end+1} = '';
        end
        outclass(unClassified) = numGroups+1;
    end
    
    % Convert back to original grouping variable
%     Edited By Mahdi
%     if isa(groupnames,'categorical')
%         labels = getlabels(groupnames);
%         if isa(groupnames,'nominal')
%             groupString = nominal(groupString,[],labels);
%         else
%             groupString = ordinal(groupString,[],getlabels(groupnames));
%         end
%         outclass = groupString(outclass);
%     elseif isnumeric(groupnames) || islogical(groupnames)
%         groupString = str2num(char(groupString)); %#ok
%         outclass = groupString(outclass);
%     elseif ischar(groupnames)
%         groupString = char(groupString);
%         outclass = groupString(outclass,:);
%     else %if iscellstr(groupnames)
%         outclass = groupString(outclass);
%     end
% End of Edition By Mahdi

else
    outclass = [];
end

outclass = - outclass;
