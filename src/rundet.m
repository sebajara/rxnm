function [t,x] = rundet(model,ts,STR1,STR1val,STR2,STR2val,STR3,STR3val,STR4,STR4val)

% Takes a MODEL struct, constructs a derivative function, then using
% its specified initial values integrates the ODE system over the
% intervals specified by TS. Return the time points T, and the species
% concentration in time X.
% 
% [T,X] = rundet(MODEL,TS)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    OPTIONAL_STRS = {'Mode','Solver','Options','Tols'};
    mode = 0; % 0 = matlab, 1 = ODEMEX
    solver = 'ode15s';
    
    options = odeset;
    [x0,p] = getinitvec(model);
    options.NonNegative = ones(size(x0));
    tols = [1e-6,1e-8,10];
    
    if nargin < 2
        error('Model must be specified');
    elseif ((nargin == 3) | (nargin == 5) | (nargin == 7) | (nargin == 9))
        error('Optional inputs must be accompanied by values');
    elseif nargin == 4
        Input_Values = {STR1, STR1val};
    elseif nargin == 6
        Input_Values = {STR1 STR1val; STR2, STR2val};
    elseif nargin == 8
        Input_Values = {STR1 STR1val; STR2, STR2val; STR3, STR3val};
    elseif nargin == 10
        Input_Values = {STR1 STR1val; STR2, STR2val; STR3, STR3val; STR4, STR4val};
    elseif nargin > 10
        error('At most 10 inputs are allowed');
    end
    for i=1:1:((nargin - 1)/2)
        STR = Input_Values{i,1};
        STRval = Input_Values{i,2};
        z = strmatch(STR,OPTIONAL_STRS,'exact');
        if z == 1
            mode = STRval;
        elseif z == 2
            solver = STRval;
        elseif z == 3
            options = STRval;
        elseif z == 4
            tols = STRval;
        else
            error ('WTH!');
        end
    end
    
    if(mode==0) % MATLAB
        [t,x] = eval([solver,'(@(t,x)',model.name,'_MATLABder(t,x,p),ts,x0,options);']);
    else
        [temp1,temp2] = eval([model.name,'_Cvode(ts,x0,p,[],tols);']);
        t = temp1';
        x = temp2';
    end

end