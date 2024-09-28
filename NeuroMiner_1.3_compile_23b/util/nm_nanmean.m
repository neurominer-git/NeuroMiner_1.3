function y = nm_nanmean(x,dim)
% FORMAT: Y = NANMEAN(X,DIM)
% 
%    Average or mean value ignoring NaNs
%
%    This function enhances the functionality of NANMEAN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEAN

% -------------------------------------------------------------------------
%    author:      Jan Gläscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $

if isempty(x)
	y = NaN;
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end

% Replace NaNs with zeros.
nans = isnan(x);
fullnans = sum(nans,dim) == size(x,2);
if anynan(x)
    try
        x(isnan(x)) = 0;
    catch ERROR
        switch ERROR.identifier
            case 'MATLAB:nomem'
                %for-loop to avoid out of memory error if x is 3-dimensional
                if length(size(x)) == 3
                    for i=1:size(x,length(size(x)))
                        temp = x(:,:,i);
                        temp(isnan(temp)) = 0;
                        try
                            x(:,:,i) = temp;
                        catch
                        end
                    end
                else
                    rethrow(ERROR)
                end
            otherwise
                rethrow(ERROR)
        end
    end
end

% denominator
count = size(x,dim) - sum(nans,dim);

% Protect against a  all NaNs in one dimension
i = find(count==0);
count(i) = ones(size(i));

y = sum(x,dim)./count;
y(i) = i + NaN;

y(fullnans) = nan;

% $Id: nm_nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
