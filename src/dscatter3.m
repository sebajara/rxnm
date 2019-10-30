function [hAxes] = dscatter3(X,Y,Z,varargin)
% DSCATTER3 creates a 3D scatter plot coloured by density.
% Based on DSCATTER by Paul H. C. Eilers and Jelle J. Goeman.
%
% [hAxes] = dscatter3(X,Y,Z,VARARGIN)
%
%   DSCATTER(X,Y,Z) creates a scatterplot of X, Y, Z at the locations
%   specified by the vectors X, Y, and Z (which must be the same
%   size), colored by the density of the points.
%    
%   DSCATTER(...,'MARKER',M) allows you to set the marker for the
%   scatter plot. Default is '+'.
% 
%   DSCATTER(...,'MSIZE',MS) allows you to set the marker size for the
%   scatter plot. Default is 10.
% 
%   DSCATTER(...,'SMOOTHING',LAMBDA) allows you to set the smoothing factor
%   used by the density estimator. The default value is 20 which roughly
%   means that the smoothing is over 20 bins around a given point.
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    lambda = [];
    nbins = [];
    plottype = 'scatter';
    %plottype = 'isosurface';
    contourFlag = false;
    msize = 10;
    marker = 'o';
    logy = false;
    filled = false;
    
    if nargin > 3
        if rem(nargin-3,2) == 1
            error('Incorrect number of arguments to %s.',mfilename);
        end
        %okargs = {'smoothing','bins','plottype','logy','marker','msize','filled'};
        okargs = {'smoothing','marker','msize'};
        for j=1:2:nargin-3
            pname = varargin{j};
            pval = varargin{j+1};
            k = strmatch(lower(pname), okargs); %#ok
            if isempty(k)
                error('Unknown parameter name: %s.',pname);
            elseif length(k)>1
                error('Ambiguous parameter name: %s.',pname);
            else
                switch(k)
                case 1  % smoothing factor
                    if isnumeric(pval)
                        lambda = pval;
                    else
                        error('Invalid smoothing parameter.');
                    end
                  case 2
                    marker = pval;
                  case 3
                    msize = pval;
                end
            end
        end
    end
    
    minx = min(X,[],1);
    maxx = max(X,[],1);
    miny = min(Y,[],1);
    maxy = max(Y,[],1);
    minz = min(Z,[],1);
    maxz = max(Z,[],1);
    
    if isempty(nbins)
        nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200), ...
                 min(numel(unique(Z)),200) ];
    end
    
    if isempty(lambda)
        lambda = 20;
    end
    
    edges1 = linspace(minx, maxx, nbins(1)+1);
    ctrs1 = edges1(1:end-1) + .5*diff(edges1);
    edges1 = [-Inf edges1(2:end-1) Inf];
    edges2 = linspace(miny, maxy, nbins(2)+1);
    ctrs2 = edges2(1:end-1) + .5*diff(edges2);
    edges2 = [-Inf edges2(2:end-1) Inf];
    edges3 = linspace(minz, maxz, nbins(3)+1);
    ctrs3 = edges3(1:end-1) + .5*diff(edges3);
    edges3 = [-Inf edges3(2:end-1) Inf];
    
    [n,~] = size(X);
    
    binyx = zeros(n,2);
    [~,binyx(:,2)] = histc(X,edges1);
    [~,binyx(:,1)] = histc(Y,edges2);
    Hyx = accumarray(binyx,1,nbins([2 1]))./n;
    Gyx = smooth1D(Hyx,nbins(2)/lambda);
    Fyx = smooth1D(Gyx',nbins(1)/lambda)';
    
    binzx = zeros(n,2);
    [~,binzx(:,2)] = histc(X,edges1);
    [~,binzx(:,1)] = histc(Z,edges3);
    Hzx = accumarray(binzx,1,nbins([3 1]))./n;
    Gzx = smooth1D(Hzx,nbins(3)/lambda);
    Fzx = smooth1D(Gzx',nbins(1)/lambda)';
    
    binzy = zeros(n,2);
    [~,binzy(:,2)] = histc(Y,edges2);
    [~,binzy(:,1)] = histc(Z,edges3);
    Hzy = accumarray(binzy,1,nbins([3 2]))./n;
    Gzy = smooth1D(Hzy,nbins(3)/lambda);
    Fzy = smooth1D(Gzy',nbins(2)/lambda)';
    
    Fyx = Fyx./max(Fyx(:));
    indyx = sub2ind(size(Fyx),binyx(:,1),binyx(:,2));
    colyx = Fyx(indyx);
    Fzx = Fzx./max(Fzx(:));
    indzx = sub2ind(size(Fzx),binzx(:,1),binzx(:,2));
    colzx = Fzx(indzx);
    Fzy = Fzy./max(Fzy(:));
    indzy = sub2ind(size(Fzy),binzy(:,1),binzy(:,2));
    colzy = Fzy(indzy);
    col = (colyx + colzx + colzy)./3;
    
    if(strcmp(plottype,'scatter'))
        siz = msize*(ones(size(col)));
        h = scatter3(X,Y,Z,siz,col,marker);
    elseif strcmp(plottype,'isosurface')
    end
    if nargout > 0
        hAxes = get(h,'parent');
    end
end

function Z = smooth1D(Y,lambda)
    [m,n] = size(Y);
    E = eye(m);
    D1 = diff(E,1);
    D2 = diff(D1,1);
    P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
    Z = (E + P) \ Y;
end