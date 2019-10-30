function [fillhandle] = myupdownfill (x,up,down,color,transparency,disc)

% This function makes a colored surface (in 2D), specified in the
% x-axis by X, and in the y-axis by UP and DOWN.  Transparency is a
% value bewteen 0 and 1, where 0 means completly transparent, and 1
% means the full color. DIST is a boolean specifying if you want to
% interpret X as discrete values.
%   
% [FILLHANDLE] = myupdownfill(X,UP,DOWN,COLOR,TRANSPARENCY,DISC)
%    
% Sebastian Jaramillo-Riveri
% June, 2013    
    
    [nr,nc] = size(x);

    if (nr > nc) % if it is a row vector, transpose everything
        x = x';
        up = up';
        down = down';
    end
    if(~(length(x)==length(up)) || (length(x)==length(down)))
        error('Dimensions not consistent');
    end
    for i = 1:length(up)
        if(~(up(i)>=down(i)))
            error('');
        end
    end

    if(disc)
        x2 = zeros(1,(length(x)*2-1));
        down2 = zeros(1,(length(x)*2-1));
        up2 = zeros(1,(length(x)*2-1));
        for i = 1:(length(x)-1)
            p = 2*i;
            x2(p:p+1) = [x(i),x(i+1)];
            down2(p:p+1) = down(i)*ones(1,2);
            up2(p:p+1) = up(i)*ones(1,2);
        end
        x2(end) = x(end);
        down2(end) = down(end);
        up2(end) = up(end);
        filled=[down2,fliplr(up2)];
        xpoints=[x2,fliplr(x2)];
    else
        filled=[down,fliplr(up)];
        xpoints=[x,fliplr(x)];
    end

    fillhandle=fill(xpoints,filled,color);
    set(fillhandle,'EdgeColor','none','FaceAlpha',transparency,'EdgeAlpha',transparency);

end