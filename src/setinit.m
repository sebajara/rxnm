function [fmodel] = setinit(imodel,id,val)
% This function takes a MODEL struct (IMODEL), and the field ID of a
% parameter or specie an set its initial condition. (Warning: it does
% not set random variables initial states).
%   
% [FMODEL] = setinit(IMODEL,ID,VAL)
%
% Sebastian Jaramillo-Riveri
% June, 2013

    fmodel = setmodel(imodel,id,'init',val);

end