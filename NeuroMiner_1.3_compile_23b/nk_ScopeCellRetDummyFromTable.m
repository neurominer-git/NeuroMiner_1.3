function [Tdummy, idxCell, idxCellInfo] = nk_ScopeCellRetDummyFromTable(act, T, uniquelim)

if ~exist("uniquelim","var") || isempty(uniquelim)
    uniquelim = 20;
end

Vars = T.Properties.VariableNames;
nVars = numel(Vars);
Tdummy = []; Vdummy = [];
idxCell = false(1,nVars);
idxCellInfo = [];       
switch act

    case 'scope'
        cnt=1;
        for i=1:nVars
            if iscell(T.(Vars{i}))
                idxCell(i)=true;
                idxCellInfo(cnt).VarName = Vars{i};
                idxCellInfo(cnt).UniqueVals = unique(rmmissing(T.(Vars{i})));
                idxCellInfo(cnt).nUniqueVals = numel(idxCellInfo(cnt).UniqueVals);
                cnt=cnt+1;
            end
        end
        
    case 'replace'
      
        for i=1:nVars
            if iscell(T.(Vars{i}))
                idxCell(i)=true;
                Ti = rmmissing(T.(Vars{i}));
                Unq = unique(Ti);
                nUnq = numel(Unq); 
                if nUnq>uniquelim
                    errordlg('\nFound %g unique values in column ''%s'' while only %g unique values are allowed.', nUnq, Vars{i}, uniquelim );
                    Tdummy = T;
                    return
                end
                Tx_dummy = nk_MakeDummyVariables(T.(Vars{i}));
                VarNames = cellstr([repmat([Vars{i} '_'], nUnq, 1) char(Unq)]);
                Tdummy = [ Tdummy array2table(Tx_dummy,"VariableNames",VarNames) ];
        
            else
                Tdummy = [ Tdummy T(:,i) ];
                Vdummy = [ Vdummy Vars{i}];
            end
        end
end