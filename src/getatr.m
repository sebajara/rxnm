function [val] = getatr(model,id,atr)

% This function takes a MODEL struct, an field string (ID) from 
% species, params, or rxns; and an atribute string (ATR) and return
% the corresponding value
%   
% [VAL] = getatr(MODEL,ID,ATR)
%   
% Sebastian Jaramillo-Riveri
% June, 2013
    
if (isfield(model.species,id))
    val = model.species.(id).(atr);
elseif(isfield(model.params,id))
    val = model.params.(id).(atr);
elseif(isfield(model.rxns,id))
    val = model.rxns.(id).(atr);
else
    error('Not a valid ID in this model');
end

end