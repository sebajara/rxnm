function [dxdt] = derivative (t,x,p)
% Derivative for the reaction model rna1
    dxdt = zeros(size(x));
    dxdt(1) = 0 -1*(p(1)*x(1)) +1*(p(2)*x(2)); % G0
    dxdt(2) = 0 +1*(p(1)*x(1)) -1*(p(2)*x(2)); % G1
    dxdt(3) = 0 +1*(p(3)*x(2)) -1*(p(4)*x(3)); % M
end
