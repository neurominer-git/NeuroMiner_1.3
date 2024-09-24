function NeuroMinerMCCMain(action, paramfile)

nm_pth = which('nm');
pth = fileparts(nm_pth);
hpc_pth = [pth filesep 'hpc'];
addpath(hpc_pth);

switch action
    case 'preproc'
        nk_Preprocess_batch(paramfile)
    case 'train'
        nk_MLOptimizer_batch(paramfile)
    case 'visualize'
        nk_VisModels_batch(paramfile)
    case 'oocv'
        nk_OOCV_batch(paramfile)
	case 'mli'
		nk_MLInterpreter_batch(paramfile)
    case 'simulate'
        nk_SimulateML_batch(paramfile)
    otherwise
        error('Action not available')
end        
