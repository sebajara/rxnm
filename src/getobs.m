function [observe] = getobs(model,obscel)

% This function takes a MODEL struct and a string cell array OBSCEL
% containing species IDs. It return a vector with the index of each
% specie in the same order.
%
% [OBSERVE] = getobs(MODEL,OBSCEL)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
observe = zeros(1,length(obscel));
for i = 1:length(obscel)
    if (isfield(model.species,obscel{i}))
        observe(i) = model.species.(obscel{i}).index;
    else
        error('String is not an species field');
    end
end

end