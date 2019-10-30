function [result] = evalinmodel(model,t,x,expr)

% Takes a MODEL struct, 
%    
% Sebastian Jaramillo-Riveri
% Aug, 2014
    [~,p] = getinitvec(model);
    %    model2eval(model);
    addpath(cd); % add path so that the function is recognized
    [result] = eval([model.name,'_eval(p,t,x,expr);']);

end