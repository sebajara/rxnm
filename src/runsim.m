function [status,result] = runsim(simo)
% This function takes a SIM structure, writes down the C code
% associated with the simulation, compiles it, and finally excecutes it. The
% messages sent to the terminal are printed. It return the final
% STATUS and RESULT comming from matlab 'system' function.
%
% [STATUS,RESULT] = runsim(simo)
%
% Sebastian Jaramillo-Riveri
% June, 2013
procsim(simo); 
[status,result] = compsim(simo);

if (status)
    result
    error('Failed compilation');
else
    if(~strcmp(result,''))
        fprintf(result);
        fprintf('\n');
    end
    excu = [simo.model.name,'_',simo.name,'_',simo.type];
    [status,result] = system(['./',excu]);
    if(~strcmp(result,''))
        fprintf(result);
        fprintf('\n');
    end
end

end