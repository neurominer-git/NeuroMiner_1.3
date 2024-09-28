function nk_PrintLogo(maindlg)
global NMinfo EXPERT SPMAVAIL DEV

clc
if ~SPMAVAIL
    %cl = [0.1,0.5,0]; 
    mode = 'non-imaging mode';
else
    %cl = NMinfo.cllogo; 
    mode = [];
end
if exist('maindlg','var') && maindlg
    fprintf('\n\t*******************************************')
    fprintf('\n\t****\\                                 /****')
    fprintf('\n\t*****\\     ');fprintf('~~~~~~~~~~~~~~~~~~~~~ '); fprintf('    /*****');
    fprintf('\n\t******\\     ');fprintf('N e u r o M i n e r ');  fprintf('    /******');
    fprintf('\n\t*******\\   ');fprintf('~~~~~~~~~~~~~~~~~~~~~ '); fprintf('  /*******');
    fprintf('\n\t*******/                           \\*******')
    fprintf('\n\t******/     ');fprintf('pattern recognition ');   fprintf('    \\******');
    fprintf('\n\t*****/      ');fprintf('for neurodiagnostic ');   fprintf('     \\*****');
    fprintf('\n\t****/           ');fprintf('applications ');      fprintf('         \\****');
    fprintf('\n\t***/                                   \\***')
    fprintf('\n\t*******************************************')
    if ~isempty(mode),fprintf('\n\t%s',mode);end
else
    fprintf('\t~~~~~~~~~~~~~~~~~~~~~~~ \n');
    fprintf('\t  N e u r o M i n e r \n');
    fprintf('\t~~~~~~~~~~~~~~~~~~~~~~~ ');
end
if EXPERT
    fprintf('\n')
    fprintf('\t>>> EXPERT MODE <<< ')
end
if DEV
    fprintf('\n')
    fprintf('\t>>> DEVELOPMENT MODE <<<')
end
fprintf('\n\t%s', NMinfo.info.ver); fprintf('\n')
if exist('maindlg','var') && maindlg
    fprintf('\n(c) %s | %s ', NMinfo.info.author, NMinfo.info.datever)
    fprintf('\n    nm@pronia.eu \n')
end