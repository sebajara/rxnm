function [fmodel] = setallspc(imodel,x0)

%   
% [FMODEL] = setallspc(IMODEL,X0)
%
% Sebastian Jaramillo-Riveri
% June, 2013

    tempmodel = imodel;

    species = fieldnames(tempmodel.species);
    nspec = length(species);
    
    for i = 1:nspec
        id = species{i};
        tempmodel = setinit(tempmodel,id,x0(i));
    end
    
    fmodel = tempmodel;
end