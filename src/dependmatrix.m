function [SDM,PDM] = dependmatrix(model)

% This function takes a MODEL structure and return two matrixes that
% relate to the dependency of each reactions rate function to
% paremeter and species values. This is infered by regular expression
% matching.
% 
% SDM is a NRXNS by NSPECIES matrix. SDM(i,j) = 0 means the
% reaction i rate does not depend on j species. 
% PDM is a NRXNS by NPARAMS matrix. PDM(i,j) = 0 means the
% reaction i rate does not depend on j parameter. 
%
% [SDM,PDM] = dependmatrix(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
% Model information
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    params = fieldnames(model.params);
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);
    
    % Instantiate output
    SDM = zeros(nrxns,nspec);
    PDM = zeros(nrxns,npara);
    
    % This regular expression represent a reparator bewteen symbols, to
    % avoid matching substrings.
    sysep = '(^|\-+|\++|\/+|\(+|\)+|\*+|\^+|$)'; % separator for symbols
                                                 % in a rxn rate function
    for i = 1:nrxns
        rxn = rxns{i};
        rate = model.rxns.(rxn).rate;
        for j = 1:nspec
            spc = species{j};
            regp = [sysep,spc,sysep];
            if (regmbool(rate,regp))
                SDM(i,j) = 1;
            end
        end
        for j = 1:npara
            par = params{j};
            regp = [sysep,par,sysep];
            if (regmbool(rate,regp))
                PDM(i,j) = 1;
            end
        end
    end
end