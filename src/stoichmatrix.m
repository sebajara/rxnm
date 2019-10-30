function [SM] = stoichmatrix(model)

% This function takes a MODEL structure and return the related
% stoichiometric matrix (SM). Row represents reactions, whereas
% columns represent species
%
% [SM] = stoichmatrix(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    nrxns = length(rxns);
    nspec = length(species);

    SM = zeros(nrxns,nspec);
    for i = 1:nrxns
        rxn = rxns{i};
        stoich = model.rxns.(rxn).stoich;
        SM(i,:) = stoich;
    end

end