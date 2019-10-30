function [fmodel] = addrvinit(imodel,syms,dist)
  
% Takes a model structure IMODEL and a string cell array SYMS, with
% species symbols. Return the same model, but seeting the initial
% value of the specie to a random variable of distribution DIST.
% DIST is formated as N by 2 matrix, where the first column are the
% values of the variable, and the second, its corresponding mass
% densety function.
%
% [FMODEL] = addrvinit(IMODEL,SYMS,DIST)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    fmodel = imodel;

    % Minor checks on DIST
    [n,m] = size(dist);
    if(~(m==2))
        error('Bad formated distribution look-up');
    end
    for i = 1:n
        if((dist(i,2)>1)||(dist(i,2)<0))
            error('Incorrect mass density values');
        end
        if(~(mod(dist(i,1),1)==0))
            error('Random initial variables should be integers');
        end
    end
    if((sum(dist(:,2))<0.9999) || (sum(dist(:,2))>1.0001))
        error('Mass densities do not add up to 1');
    end
    for i = 1:length(syms)
        sym = syms{i};
        if (isfield(fmodel.species,sym))
            fmodel.species.(sym).type = 'RV2';
            fmodel.species.(sym).dist = dist;
            fmodel.species.(sym).init = distmoments(dist);
        elseif(isfield(fmodel.params,sym))
            fmodel.params.(sym).type = 'RV2';
            fmodel.params.(sym).dist = dist;
            fmodel.params.(sym).init = distmoments(dist);
        else
            error('Not a valid symbol in the model');
        end
    end

end