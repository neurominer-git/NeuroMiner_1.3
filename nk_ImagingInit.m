function [spmrootdir, fsrootdir, jurootdir] = nk_ImagingInit(neurominerpath, imaging_init_path, delete_existing)
global DEV
if ~exist('delete_existing','var') || isempty(delete_existing)
    delete_existing = false;
end

if ~exist(imaging_init_path,'file') || delete_existing
    spmrootdir = uigetdir(neurominerpath,'Specify SPM directory');
    fsrootdir = uigetdir(neurominerpath,'Specify Freesurfer MATLAB directory');
    if DEV
        jurootdir = uigetdir(neurominerpath,'Specify JuSpace directory');
    else
        jurootdir = [];
    end
    save(imaging_init_path,'spmrootdir','fsrootdir','jurootdir');
else
    P = load(imaging_init_path); 
    spmrootdir=P.spmrootdir;
    fsrootdir=P.fsrootdir;
    jurootdir=[];
    if isfield(P,'jurootdir'), jurootdir=P.jurootdir; end
    % Check whether paths exist and update paths accordingly
    if exist('spmrootdir','var') && ~isnumeric(spmrootdir) && ~exist(spmrootdir,'dir') 
        tspmrootdir = fileparts(which('spm'));
        if ~isempty(tspmrootdir) 
            spmrootdir = tspmrootdir; 
            try
                save(imaging_init_path,'spmrootdir','-append')
                fprintf('\nUpdated the SPM paths according to the system settings.')
            end
        end
    end
    if exist('fsrootdir','var') && ~isnumeric(fsrootdir) && ~exist(fsrootdir,'dir')
        tfsrootdir = fileparts(which('MRIread'));
        if ~isempty(tfsrootdir) 
            fsrootdir = tfsrootdir; 
            try
                save(imaging_init_path,'fsrootdir','-append')
                fprintf('\nUpdated the SPM paths according to the system settings.')
            end
        end
    end
    if DEV
        if exist('jurootdir','var') && ~isnumeric(jurootdir) && ~exist(jurootdir,'dir')
            tjurootdir = fileparts(which('JuSpace'));
            if ~isempty(tjurootdir)
                jurootdir = tjurootdir;
                try
                    save(imaging_init_path,'jurootdir','-append')
                    fprintf('\nUpdated the JuSpace toolbox path according to the system settings.')
                end
            end
        end
    end
end