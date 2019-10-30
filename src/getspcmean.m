function [xs,stds] = getspcmean(model)
%
%
% [XS,STDS] = getspcmean(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    species = fieldnames(model.species);
    nspec = length(species);

    xs = zeros(1,nspec);
    stds = zeros(1,nspec);
    for i = 1:nspec
        id = species{i};
        typ = model.species.(id).type;
        ind = model.species.(id).index;
        if (strcmp(typ,'RV1'))
            xs(ind) = model.species.(id).init;
        elseif((strcmp(typ,'RV2') || strcmp(typ,'RV3')))
            [x1,x2] = distmoments(model.species.(id).dist);
            xs(ind) = x1;
            stds(ind) = sqrt(x2-x1^2);
        end
    end

end