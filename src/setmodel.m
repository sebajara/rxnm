function [fmodel] = setmodel(imodel,id,atr,val)

% This function takes a MODEL struct (IMODEL), and the field ID of a
% parameter,specie or rxn, an set its atribute (ATR) to VAL.
%   
% [FMODEL] = setmodel(IMODEL,ID,ATR,VAL)
%
% Warning! It does not perform any special type checking
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
tempmodel = imodel;
if (isfield(tempmodel.species,id))
    tempmodel.species.(id).(atr) = val;
elseif(isfield(tempmodel.params,id))
    tempmodel.params.(id).(atr) = val;
elseif(isfield(tempmodel.rxns,id))
    tempmodel.rxns.(id).(atr) = val;
else
    error('Not a valid specie or parameter ID in this model');
end
fmodel = tempmodel;

end