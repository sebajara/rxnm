function [caf] = sortcabylength(cai,order)
% 
% Sort a string cell array (CAI) by the lenght of the string
% 
% [CAF] = sortcabylength(CAI)
% 
% Sebastian Jaramillo-Riveri. 
% June, 2013
    
    if(iscellstr(cai))
        n = length(cai); % cell array size
        caf = cell(n,1);
        lengths = zeros(n,1);
        for i = 1:n
            lengths(i) = length(cai{i}); 
        end
        if(order)
            [~,sls] = sort(lengths,'ascend');
        else
            [~,sls] = sort(lengths,'descend');
        end
        caf = cai(sls);
    else
        error('Argument is not a string cell array\n');
    end

end