function [Dummy, Dummynum, Fu] = nk_MakeDummyVariables(V,vec,mode)

if ~exist('mode','var'), mode = 'stable'; end
[m, n] = size(V);
if m > 1 && n > 1, error('Only vector operations are supported!'); end
if exist('vec','var') && ~isempty(vec)
    Fu = vec; nFu = numel(vec);
else
    Fu = unique(rmmissing(V),mode); nFu = numel(Fu);
end
if m > 1
    Dummy = false(m, nFu);
else
    Dummy = false(nFu, n);
end

for i = 1:nFu
    try
        ind = V == Fu(i);
    catch
        ind = strcmp(V,Fu{i});
    end
    if m > 1 
        Dummy(ind,i) = true;
    else
        Dummy(i, ind) = true;
    end
end
[~, Dummynum] = max(Dummy,[],2); 
Ix = sum(Dummy,2)==0;
if any(Ix)
    Dummy = double(Dummy);
    Dummy(Ix,:) = NaN;
    Dummynum(Ix)=NaN;
end

