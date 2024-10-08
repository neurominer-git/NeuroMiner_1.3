function param = nk_SVM_config(res, param, progtype, kerntype, parentstr)

if ~exist('progtype','var'),   progtype = []; end
if ~exist('kerntype','var'),   kerntype = []; end

%% Define SVM optimization options depending on the SVM implementation
if isempty(progtype)
     sftmenu = ['SVM --------------> LIBSVM|' ...
                'SVM --------------> LIBLINEAR'];
            
     sftval = ['LIBSVM';'LIBLIN'];
     progtype = nk_input(['Select prediction method for ' res.modeflag ' framework'],0,'m', sftmenu, sftval, 1);
     param.SVM.prog = progtype;
end

switch progtype
    
    case 'LIBSVM'
        
        if isfield(res.TrainParam,'LABEL') && res.TrainParam.LABEL.flag == 1
            modeflag = res.TrainParam.LABEL.newmode;
        else
            modeflag = res.modeflag;
        end
        act = 1; 
        while act
            [act,param] = nk_LIBSVM_config(res, param, [], [], parentstr, modeflag); 
        end
        
    case 'SVMPRF'
        
        param = nk_SVMPRF_config(param);
        
    case 'LIBLIN'
        
        if isfield(res.TrainParam,'LABEL') && res.TrainParam.LABEL.flag == 1
            modeflag = res.TrainParam.LABEL.newmode;
        else
            modeflag = res.modeflag;
        end
        act = 1;
        while act
            [act, param] = nk_LIBLIN_config(modeflag, param, [], parentstr);
        end
        kerntype = [];
        
    case 'LSTSVM'
        
        param = nk_LSTSVM_config(param);
            
end
%% Define SVM kernel
if isempty(kerntype)
    param = nk_Kernel_config(param, 1);
else
    param.kernel.kernstr = kerntype;
end        

end