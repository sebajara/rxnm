function [figs] = ploscatter(simo,scatters)
% Takes a 'ss' type SIM structure, and assuming it was completed it
% will parse the output and generate density scatter plots specified
% pair or triplet of species IDs. Each element in the cell array
% SCATTERS is a string cell array of 2 or 3 elements containing
% species IDs.  Figures in .fig format will be saved.
%     
% [FIGS] = ploscatter(simo,SCATTERS)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
name = simo.name;
model = simo.model;
type = simo.type;
ts = simo.ts;
nrep = simo.nrep;

if(~strcmp(type,'ss'))
    error('plotscatter is meant to be used for ''ss'' type simulation outputs');
end

xtitles = model.ind2spc;
for i = 1:length(xtitles)
    xtitles(i) = regexprep(xtitles(i),'_','\\_');
end

[nscatter,~] = size(scatters);
ntps = length(ts);

[T,X] = loadsimout(simo,model.ind2spc); % load all data once

for i = 1:nscatter
    scatter = scatters{i};
    observe = getobs(model,scatter);
    nobs = length(observe);
    if((nobs>3)||(nobs<2))
        error('Invalid number of species for scatter plots');
    end
    fig = figure(); hold on;
    if(nobs==2)
        j1 = observe(1) + 1;
        j2 = observe(2) + 1;
        xs1 = X{j1-1,1};
        xs2 = X{j2-1,1};
        max1 = 0;
        max2 = 0;
        ms1 = zeros(ntps,1);
        ms2 = zeros(ntps,1);
        for tp = 1:ntps
            subplot(ntps,1,tp); hold on;
            x1 = xs1(tp,:);
            x2 = xs2(tp,:);
            if((std(x1)==0)||(std(x2)==0))
                plot(x1,x2,'+');
            else
                %plot(x1,x2,'+','LineWidth',2);
                dscatter(x1',x2','MARKER','+');
                %dscatter(x1',x2','MSIZE',20);
                %dscatter(x1',x2','plottype','contour');
            end
            ms1(tp) = mean(x1);
            ms2(tp) = mean(x2);
            if(max(x1)>max1)
                max1 = max(x1);
            end
            if(max(x2)>max2)
                max2 = max(x2);
            end
        end
        for tp = 1:ntps
            subplot(ntps,1,tp); hold on;
            plot([0 max1],[ms2(tp) ms2(tp)],':','Color','black','LineWidth',2);
            plot([ms1(tp) ms1(tp)],[0 max2],':','Color','black','LineWidth',2);
            set(gca,'FontSize',20);
            xlabel(xtitles{j1-1},'FontSize',25); ylabel(xtitles{j2-1},'FontSize',25);
            title(['t=',num2str(T(tp))],'FontSize',25);
            axis([0 (max1+1) 0 (max2+1)]);
        end
        saveas(fig,[model.name,'_',name,'_ss',xtitles{j1-1},'_vs_',xtitles{j2-1}], 'fig');
        figs{i} = fig;
    elseif(nobs==3)
        j1 = observe(1) + 1;
        j2 = observe(2) + 1;
        j3 = observe(3) + 1;
        xs1 = X{j1-1,1};
        xs2 = X{j2-1,1};
        xs3 = X{j3-1,1};
        max1 = 0;
        max2 = 0;
        max3 = 0;
        ms1 = zeros(ntps,1);
        ms2 = zeros(ntps,1);
        ms3 = zeros(ntps,1);
        for tp = 1:ntps
            subplot(ntps,1,tp); hold on;
            x1 = xs1(tp,:);
            x2 = xs2(tp,:);
            x3 = xs3(tp,:);
            if((std(x1)==0)||(std(x2)==0)||(std(x3)==0))
                plot3(x1,x2,x3,'+');
            else
                ms1(tp) = mean(x1);
                ms2(tp) = mean(x2);
                ms3(tp) = mean(x2);
                dscatter3(x1',x2',x3');
                %dscatter(x1',x2');
                %dscatter(x1',x2','PLOTTYPE','contour');
            end
            if(max(x1)>max1)
                max1 = max(x1);
            end
            if(max(x2)>max2)
                max2 = max(x2);
            end
            if(max(x3)>max3)
                max3 = max(x3);
            end
        end
        for tp = 1:ntps
            subplot(ntps,1,tp); hold on;
            %view(45,45);
            view(-60,60)
            set(gca,'FontSize',15);
            xlabel(xtitles{j1-1},'FontSize',20); ylabel(xtitles{j2-1},'FontSize',20);
            zlabel(xtitles{j3-1},'FontSize',20);
            %title(['t=',num2str(T(tp))],'FontSize',20);
            axis([0 (max1+1) 0 (max2+1) 0 (max3+1)]);
            grid on;
        end
        saveas(fig,[model.name,'_',name,'_ss',...
                    xtitles{j1-1},'_vs_',xtitles{j2-1},...
                    '_vs_',xtitles{j3-1}], 'fig');
        figs{i} = fig;
    end
end

end