function [sY, IN] = cv_perfROImeans(Y, IN)%, PREPROC, prevP)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfROImeans(Y{i}, IN); end %, PREPROC, prevP); end
    else
        [ sY, IN ] = PerfROImeans(Y, IN);%, PREPROC, prevP);
    end
end

% =========================================================================
function [Y, IN] = PerfROImeans(Y, IN)%, PREPROC, prevP)
    
    inputStr = '';
    
    if isfield(IN, 'V_brainmask') && isfield(IN, 'atlas') && isfield(IN, 'YAtlas') && isfield(IN, 'indVol_brainmask')
        [Y, IN] = cv_compute_ROImeans(Y, IN);%, PREPROC, prevP);
    end
    
end