function [CA] = mysplit(str,sep)
% Function to splist a string STR using a regular expression (SEP).
% For some reason the split version of my MATLAB does not support
% regex. Returns up to 100 values (arbitrary)
% 
% [CA] = mysplit(STR,SEP)
%
% Sebastian Jaramillo-Riveri
% June, 2013        
    
    maxc = 100;
    nmatch = 0;
    mstr = ['(.+?)','(?:^|',sep,'|$)','(.*)'];
    CA = cell(maxc,1);
    
    r = str;
    while (regmbool(r,mstr) && (nmatch<=maxc))
        t = regexp(r,mstr,'tokens'); % Find the first match
        t = t{1}; % Adjust the cell to the ouput of regexp tokens
        if (~regmbool(t{1},sep)) % Discard matches to the separator
            nmatch = nmatch + 1;
            CA{nmatch} = t{1}; % Save the first match
        end
        if (length(t)>1)
            r = t{2}; % Use the rest for the next iteration
        else
            r = '';
        end
    end

    if (nmatch)
        CA = CA(1:nmatch);
    else
        CA = {};
    end

end