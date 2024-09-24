function [NM, Y, oocvind, fldnam, dattype] = nk_SelectOOCVdata(NM, oocvflag, selflag, multiflag)

oocvind = []; fldnam = []; dattype = [];
if ~exist('selflag','var') || isempty(selflag), selflag = 0; end
if ~exist('multiflag','var') || isempty(multiflag), multiflag = 0; end
if ~isfield(NM,'C') && ~isfield(NM,'OOCV') && ~selflag, return; end
switch selflag 
    case 0
        selstr = 'select';
    case 1
        selstr = 'modify';
end
Y = [];
switch oocvflag
    case 1
        fldnam = 'OOCV'; St = 'Independent Test Data Manager'; dattype = 'independent test data';
    case 2
        fldnam = 'C'; St = 'Calibration Data Manager'; dattype = 'calibration data';
end
if isfield(NM,fldnam), Y = NM.(fldnam); end 
act = true;

while act  

    O = nk_OOCVDataIO('title',St,'list',Y,'mode',selstr,'multiselection',multiflag);
    
    if isfield(O,'delete')
        Y(O.delete)=[];
        NM.(fldnam)(O.delete)=[];
    elseif isfield(O,'Items') && ~isempty(O.Items)
        NM.(fldnam) = O.Items;
    end
    if isfield(O,'exit') && O.exit
        act = false;
    end
end
if isfield(O,'SelItem')
    if ~isempty(O.SelItem)
        oocvind = O.SelItem;
    else
        oocvind = 1; 
    end
    Y = O.Items{oocvind};
end
