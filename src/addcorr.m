function [fmodel] = addcorr(imodel,syms,c)

% Takes a model structure IMODEL and a string cell array SYMS, with
% species symbols. Return the same model, but seeting the correlation
% bewteen all species to C.
%
% [FMODEL] = addcorr(IMODEL,SYMS,C)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
fmodel = imodel;

if(length(syms)<2)
    error('At least two simbols are needed to set correlations');
end
if((c<0)||(c>1))
    error('Correlations should be bewteen 0 and 1');
end


for i = 1:length(syms)
    sym1 = syms{i};
    for j = 1:length(syms)
        sym2 = syms{j};
        if (~strcmp(sym1,sym2)) % not equal
            if (isfield(fmodel.species,sym1) && ...
                isfield(fmodel.species,sym2) )
                if(isfield(fmodel.species.(sym1),'dist') && ...
                   isfield(fmodel.species.(sym2),'dist') )
                    ind1 = fmodel.species.(sym1).index;
                    ind2 = fmodel.species.(sym2).index;
                    fmodel.corr(ind1,ind2) = c;
                    fmodel.corr(ind2,ind1) = c;
                else
                    error('Species dont have a distribution');
                end
            else
                error('Not a valid symbol species in the model');
            end
        end
    end
end

end