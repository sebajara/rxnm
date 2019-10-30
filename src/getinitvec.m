function [x0,p] = getinitvec(model)
% Takes a MODEL struct and returns the initial state values for each
% specie (X0) and the value of parameters (P).
%
% [X0,P] = getinitvec(MODEL)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
species = fieldnames(model.species);
params = fieldnames(model.params);

nspec = length(species);
npara = length(params);

x0 = zeros(1,nspec);
p = zeros(1,npara);

for i = 1:nspec
    id = species{i};
    val = model.species.(id).init;
    x0(i) = val;
end
for i = 1:npara
    id = params{i};
    val = model.params.(id).init;
    p(i) = val;
end

end