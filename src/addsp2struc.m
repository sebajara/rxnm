function [fmodel,nf] = addsp2struc(imodel,kind,ni,sym,val)
% Parameter and species have basically the same structure. This
% function initialize them, and makes sure they are unique without
% replicates.
% imodel := model structure
% kind := either 'species' or 'params' (string)
% ni := current number of that particular type. Index would be that
%       value plus 1
% sym := name 
% val := init value
    fmodel = imodel;
    nf = ni;
    if (~isfield(fmodel.sym2ind,sym))
        nf = nf + 1;
        fmodel.('sym2ind').(sym) = nf;
        fmodel.(kind).(sym) = struct;
        fmodel.(kind).(sym).name = sym;
        fmodel.(kind).(sym).index = nf;
        if(strcmp(kind,'vars'))
            fmodel.(kind).(sym).formula = val;
        else
            fmodel.(kind).(sym).init = val;
        end
        fmodel.(kind).(sym).('type') = 'RV1'; % for stochastic interpretation only
    end
end
