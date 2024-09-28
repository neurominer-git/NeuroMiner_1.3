function handles = load_analysis(handles, varargin)

% Loop through input params
nVarIn = length(varargin);

for i = 1:nVarIn

    if strcmpi(varargin{i}, 'Subjects')
        
        handles.subjects = varargin{i+1};
    
    elseif strcmpi(varargin{i}, 'Params')
        
        handles.params = varargin{i+1};
    
    elseif strcmpi(varargin{i}, 'Analysis')
        
        handles.GDdims = varargin{i+3};
        handles = load_GDdims(handles, varargin{i+1}, handles.label, handles.GDdims);

    elseif strcmpi(varargin{i}, 'Visdata')
    
        vis = varargin{i+1};
        if ~isempty(vis)
            [nM, nL] = size(vis);
            for m=1:nM
                for l=1:nL
                    handles.visdata_table(m, l) = create_visdata_tables(vis{m,l}, [], [], 'create');
                end
            end
            handles.visdata = vis;
        end
        
    elseif strcmpi(varargin{i}, 'OOCVdata')
        
        handles = load_OOCV(handles, varargin{i+1});

    elseif strcmpi(varargin{i}, 'MLIdata')
        
        mli = varargin{i+1};
        handles.MLIdata = mli; 
        if ~isfield(handles, 'MLIapp') 
            handles.MLIapp = 0; 
        end
    end
    
    
end
