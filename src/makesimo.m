function [simo] = makesimo(model,name,type,ts,nrep)
% Shortcut to instantiate a SIM structure
%
% [SIMO] = makesimo(MODEL,NAME,TYPE,TS,NREP)
% 
% Sebastian Jaramillo-Riveri
% June, 2013

simo = struct;
simo.model = model;
simo.type = type;
simo.ts = ts;
simo.nrep = nrep;
simo.name = name;

end