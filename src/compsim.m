function [status,result] = compsim(simo)
    
% This function takes a SIM structure, and compiles the corresponding
% C code, given that it exists! It follows the naming conventions from
% the function 'procjob'. It return the STATUS and RESULT from matlab
% 'system' function.
% 
% [STATUS,RESULT] = compsim(SIMO)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    name = simo.name;
    model = simo.model;
    type = simo.type;
    ts = simo.ts;
    nrep = simo.nrep;
    
    compiler = 'gcc'; % GNU C compiler
    %extra = 'libmeschach.a';
    extra = '';
    % Naming convention
    if (strcmp(type,'ss'))
        cfile = [model.name,'_',name,'_ss.c'];
        ofile = [model.name,'_',name,'_ss'];
    elseif(strcmp(type,'st'))
        cfile = [model.name,'_',name,'_st.c'];
        ofile = [model.name,'_',name,'_st'];
    end
    
    tocomp = [compiler,' ',cfile,' ',extra,' -o ',ofile,' -lm'];
    
    [status,result] = system(tocomp);

end