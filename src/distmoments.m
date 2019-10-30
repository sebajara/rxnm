function [x1,x2] = distmoments(dist)

% This function takes a NVALUES by 2 matrix (DIST), and interprets it
% a mass densety function, where the first columns represent the
% values of the random variable, and the second colum the
% corresponding mass density. It returns the first and second
% moments (X1, and X2).
%
% [X1,X2] = distmoments(DIST)
%
% Sebastian Jaramillo-Riveri
% June, 2013

x1 = sum(dist(:,1).*dist(:,2));
x2 = sum((dist(:,1).^2).*dist(:,2));

end