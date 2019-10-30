function [combs,pointind] = combsets(CA)

% Takes a cell array containing vectors and returns a matrix with all
% possible combinations bewteen their elements (COMBS). Where each row
% is a combination, and columns correspond to order of vectors in the
% cell array. Also returns a matrix with the indexes of each value
% in COMBS to the arrays in CA
%
% [COMBS] = combsets(CA)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
[n,m] = size(CA);
if(n>m)
    nsets = n;
else
    nsets = m;
end

nvals = zeros(nsets,1);
for i = 1:nsets
    nvals(i) = length(CA{i});
end

ncombs = prod(nvals(1:nsets));
combs = zeros(ncombs,nsets);
pointind = zeros(ncombs,nsets);

for s = 1:nsets
    nrep = prod(nvals(s:nsets))/nvals(s);
    t = ncombs/(nrep*nvals(s));
    i = 1;
    vals = CA{s};
    for k = 1:t
        for j = 1:nvals(s)
            combs(i:(i+nrep-1),s) = vals(j).*ones(nrep,1);
            pointind(i:(i+nrep-1),s) = j.*ones(nrep,1);
            i = i + nrep;
        end
    end
end

end