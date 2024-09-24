function  [NM, act] = nk_DataIO_config( NM, parentstr, act, varind, IO )

mn_str=[];
nk_PrintLogo
if isfield(NM,'Y') && ~isempty(NM.Y), Yfl = 1; else, Yfl = 0; end
if isfield(NM,'covars') && ~isempty(NM.covars), covfl = 1; else, covfl = 0; end
mestr = 'NM data workspace manager'; navistr = [parentstr ' >>> ' mestr];

if ~exist('act','var') || isempty(act)
    if ~Yfl
        act = 1;
    else
        nk_SelectVariateIndex(NM);
        mn_str = [ mn_str 'Add modality to NM workspace|' ...
                          'Modify modality in NM workspace|' ...
                          'Delete modality from NM workspace|' ...
                          'Delete all modalities from NM workspace|' ...
                          'Open NM Modality Manager ...']; 
        mn_act = [ 2 3 4 5 9 ];
        if ~covfl 
            fprintf('\n\n');
            fprintf('--- No covariates found in NM workspace.');
            mn_str = [mn_str '|Add covariate(s) to NM workspace' ];
            mn_act = [mn_act 6];
        else
            nk_SelectCovariateIndex(NM, [], 0);
            mn_str = [mn_str '|Modify covariate(s) in NM workspace|Delete covars from NM workspace|' ];
            mn_act = [mn_act 6 7];
        end
        mn_str = [mn_str '|Finish data import for discovery & cross-validation analysis|' ]; mn_act = [mn_act 8];
        fprintf('\n\n'); fprintf('You are here: %s >>> ',parentstr); 
        act = nk_input(mestr,0,'mq', mn_str, mn_act);
    end
end

switch act 

    case 0
        return
        
    case {1,2,3} % Setup NM workspace for analysis
        switch act
            case 1
                 % Setup NM workspace and import data to new modality
                 varind = 1; IO=[];
                 if isfield(NM,'datadescriptor') && isfield(NM.datadescriptor{varind},'input_settings')
                     IO = NM.datadescriptor{varind}.input_settings;
                 end
            case {2,3}
                 % Import data to new NM modality
                if ~exist('varind','var') || isempty(varind)
                    switch act 
                        case 2
                            varind = numel(NM.Y)+1;
                        case 3
                            varind = nk_SelectVariateIndex(NM, 1, [], 1);
                    end
                end
                if  ~exist('IO','var') || isempty(IO)
                    if varind <= numel(NM.datadescriptor)
                        IO = NM.datadescriptor{varind}.input_settings;
                    else
                        IO = [];
                    end
                end
        end
        t_act = 'loop'; t_mess = [];while ~strcmp(t_act,'BACK'), [ IO, t_act, t_mess ] = DataIO( NM, mestr, IO, t_mess, varind) ; end
        if IO.completed
            [ NM, NM.datadescriptor{varind}.input_settings ] = TransferModality2NM( NM, IO, varind ); 
        else
            NM.datadescriptor{varind}.input_settings = IO;
        end
        if isfield(NM,'Y') && ~isempty(NM.Y), Yfl = true; else, Yfl = false; end
        if ~Yfl, act = 0; end
        
    case 4
        varind = nk_SelectVariateIndex(NM, 1, [], 1);
        delfl = nk_input(sprintf('Do you really want to delete modality %g from the NM workspace',varind),0,'yes|no',[1,0],2);
        if delfl
            NM = DeleteModalityInNM(NM, varind);
            if isempty(NM.datadescriptor)
                NM = DeleteAll(NM);
            end
        end

    case 5
        delfl = nk_input('Do you really want to completely reset the NM workspace',0,'yes|no',[1,0],2);
        if delfl
             NM = DeleteAll(NM);
        end

    case 6
        if covfl
            t_covars = NM.covars;
        else
            t_covars = [];
        end
        [NM.covars, NM.covnames] = nk_DefineCovars_config(NM.n_subjects_all, t_covars);
         
    case 7
        covind = nk_SelectCovariateIndex(NM,1,1);
        delfl = nk_input(sprintf('Do you really want to delete the covariates from the NM workspace'),0,'yes|no',[1,0],2);
        if delfl, NM = DeleteCovarsInNM(NM, covind); end
        
    case 8
        NM.id = genID;
        NM.defs.import_finished = true; act = 0;

    case 9
        inp.currmodal = 1;
        inp.nummodal = numel(NM.Y);
        inp.na_str = '?';
        inp.covflag = false;
        inp.desc = NM.datadescriptor{inp.currmodal}.desc;
        inp.Ydims = cellfun('size', NM.Y, 2);
        inp.datadescriptor = NM.datadescriptor;
        inp.brainmask = NM.brainmask;
        inp.badcoords = NM.badcoords;
        inp.ncases = size(NM.cases,1);
        inp.oocvmode = false;
        if isfield(NM,'id')
            inp.id = NM.id;
        end
        inp.modeflag = NM.modeflag;
        if isfield(NM,'covars'), inp.covflag = true; inp.covars = NM.covars; end

        actx = 1; 
        while actx
            [actx, NM, inp] = nk_ManageDataContainers(actx, inp, NM, parentstr); 
        end

end

