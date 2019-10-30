function [fig] = plotsim(simo,obscel)

% This function takes a SIM structure and a string cell array
% containing species names (OBSCEL). Assuming that the simulation was
% 'executed', it will parse the output files and obtain the
% corresponding data for the species specified in OBSCEL, and plot the
% results. It will save the figure in .fig format
%
% For 'ss' simulations, it will plot the cumulative frequency of each
% observed specie. Each color will represent a time point.
%
% For 'st' simulations, it will plot the average value of each specie.
% Each color will represent one simulation run.
%     
% [fig] = plotssim(SIMO,OBSCEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
name = simo.name;
model = simo.model;
type = simo.type;
ts = simo.ts;
nrep = simo.nrep;

xtitles = model.ind2spc;
for i = 1:length(xtitles)
    xtitles(i) = regexprep(xtitles(i),'_','\\_');
end
observe = getobs(model,obscel);
nobs = length(observe);
ntps = length(ts);
transparency = 0.5;

if (strcmp(type,'ss'))
    colors = hsv(ntps);
    [T,X] = loadsimout(simo,obscel);
    fig = figure(); hold on;
    for obs = 1:nobs
        maxX = 0;
        o = observe(obs);
        j = o + 1;
        subplot(nobs,1,obs); hold on;
        xs = X{obs,1};
        for tp = 1:ntps
            x = xs(tp,:);
            h = myrelfreq(x,[0:max(x)]);
            stairs([0:max(x)],cumsum(h),'Color',colors(tp,:),'LineWidth',2);
        end
        if (max(x) > maxX)
            maxX = max(x);
        end
        set(gca,'FontSize',20);
        xlabel(xtitles{j-1},'FontSize',25); ylabel('c.m.f','FontSize',25);
        axis([0 (maxX+1) 0 1]);
    end
    saveas(fig,[model.name,'_',name,'_ss'], 'fig');
elseif(strcmp(type,'st'))
    colors = hsv(nrep);
    [T,X,STD] = loadsimout(simo,obscel);
    fig = figure(); hold on;
    t = T(:,1) + (T(:,2) - T(:,1))./2;
    for obs = 1:nobs
        maxX = 0;
        o = observe(obs);
        j = o + 2;
        subplot(nobs,1,obs); hold on;
        xs = X{obs,1};
        stds = STD{obs,1};
        for r = 1:nrep
            x = xs(:,r);
            std = stds(:,r);
            if (max(x)>maxX)
                maxX = max(x);
            end
            stairs(t,x,'Color',colors(r,:),'LineWidth',2);
            %myupdownfill(t,x+2*std,x-2*std,colors(r,:),transparency,1);
        end
        set(gca,'FontSize',12);
        axis([0 ts(end) 0 (maxX+1)]);
        xlabel('time','FontSize',20); ylabel(xtitles{o},'FontSize',20);
    end
    saveas(fig,[model.name,'_',name,'_st'], 'fig');
end

end