function [orate] = odemexsynrate(rate)

% Take a rate expression, and fix the syntax such that is compatible
% for compilation with ODEMEX
% -> Change ^ operator by 'pow'

% Work only with parenthesis
    
    if(regmbool(rate,'\^'))
        sides = strsplit(rate,'^');
        side1 = sides{1};
        side2 = sides{2};
        if(strcmp(side1(end),')'))
            parcount = 0;
            parend = 0;
            nside1 = length(side1);
            for j = 1:nside1-1
                ch = side1(nside1-j);
                if(strcmp(ch,')'))
                    parcount = parcount-1;
                elseif(strcmp(ch,'('))
                    parcount = parcount+1;
                    if(parcount==1) % right one
                        parend = nside1-j;
                        break;
                    end
                end
            end
        end
        %% will not work everytime
        orate = [side1(1:parend-1),'pow(',side1(parend:end),',',side2,')'];
    end
end