function [] = parsejob(file)
% Parse a file with .job extension and process it. After the
% simulations are done, it will create plots and saves the data and
% figures to a file. The prefix for most files is 
% <MODEL.name>_<NAME>_
% All lines beginning with # are interpreted as comments.
% 
% Valid options are:
%
% NAME := will be used as part of the files created later. Optional.
% Example:
% NAME test1
%
% MODEL := file in .rxnm or .sbproj with the reaction model
% definition. Mandatory. 
% Example:
% MODEL test.rxnm
%
% OUTFOLDER := folder where the subsequent files will be created If
% the folder does not exists it will be created. Optional. 
% Example: 
% OUTFOLDER simtest1/
%
% DOT := If is 0, it won't attempt to create a graph representation of
% the model. Requires graphviz to be installed. Optional. 
% Example:
% DOT 1
%
% SETINIT := symbol1,symbol2,...,symboln:value. It will set the
% initial value of every symbol (parameter or specie) to value. Also
% symbols can be specified by regular expressions. Optional. 
% Example:
% SETINIT k1:10
% SETINIT x1,x2,x3:10
% SETINIT /^x\d+$/:10
% 
% SETDIST := symbol1,symbol2,...,symboln:file.dat. For each specie
% specified by symbol, it sets its initial condition to a random
% variable with mass density function defined in file.dat. file.dat is
% a text file with 2 columns, the first are species state and the
% second the corresponding mass density. Optional.
% Example:
% SETDIST x1,x2,x3:poiss10.dat
% SETDIST /^x\d+$/:poiss10.dat
%
% SETCORR := symbol1,symbol2,...,symboln:value. It will set the
% correlation beteen all species specified by symbols to
% value. value should be between 0 and 1, and species should have
% an associated distribution for their initial condition.
% Example:
% SETCORR x1,x2,x3:0.8
% SETCORR /^x\d+$/:0.8
%
% TYPE := Simulation jobs. They can be separated by commas. Valid
% values are: DET,SS,ST,SCANDET,SCANSS.
% DET represents deterministic simulation, it plots the trayectories
% for each specie in OBSERVE and saves the data and figure (_det at
% the end).
% SS performs stochastic simulations computing state distributions,
% plotting the cumulative frequency of each specie in OBSERVE at each
% time point, and saving those series and figures (_ss at the end).
% ST performs stochastic simulations computing the average species
% states in time intervals, for each simulation it plots the
% trayectorie of the mean over time, saves the series and figures (_st
% at the end).
% SCANDET scans the values of species at time SCANTIME, variating
% parameters and/or initial states as specified by VARIATE. Output
% will be saved with the extension _scandet.
% SCANSS scans the states of species time SCANTIME, variating
% parameters and/or initial states as specified by VARIATE, taking
% SCANREP samples per point.
% Optionals.
% Example:
% TYPE DET,ST
% TYPE SS
% 
% OBSERVE := symbol1,symbol2,...,symboln. After the simulations are
% done, we will analyse and plot only the species specified here. When
% not specified, it will run the simulations, but will not make
% any plots.
% Examples: 
% OBSERVE x1,x2,x3
% OBSERVE ALL
%
% SCATTER := A couple or a triplet of species IDs separated by
% commas. A scatter plot of the states of those species at each
% point from SSTPOINTS would be made, and saved in .fig
% format. 
% Examples:
% SCATTER x1,x2
% SCATTER x1,x2,x3
%
% DETPOINTS := time points used for numerical integration in the
% deterministic setting. Syntax is the same as for defining vectors
% in MATLAB. Required for DET job only. 
% Example:
% DETPOINTS 0:1:1000 
%
% SSTPOINTS := time points for which the state of the system will be
% saved during simulations. Syntax is the same as for defining vectors
% in MATLAB. Once plotted, each time point will be represented by a
% color. Required for SS job only. 
% Example:
% SSTPOINTS 10,100,1000
%
% STTPOINTS := time intervals for which the state of the species
% will be averaged. Syntax is the same as for defining vectors
% in MATLAB. Required for ST job only. 
% Example:
% STTPOINTS 0:1:1000
%
% SSNREP := Number of times we should simulate in the SS type of job,
% to compute states distribution. Required for SS job only. 
% Example:
% SSNREP 1000
% 
% STNREP := Number of times we should simulate in the SS type of
% job. Once plotted, each trajectory will be represented by a
% color. Required for ST job only. 
% Example:
% STNREP 10
% 
% SCANTIME := time point for SCANDET and SCANSS. Required for both
% of those job types.
% Example:
% SCANTIME 100000
%
% VARIATE := specifies the points of parameters or initial states for
% SCANDET and SCANSS. Won't have any effect on initial states sampled
% from distributions. When it appears more than once, it will scan
% every combination of values including the preceding ones. Required
% for SCANDET and SCANSS.
% Example:
% VARIATE x 0:10
% VARIATE k 0:0.1:1
% 
% SCANREP := Specifies the number of samples for the stochastic
% simulations in SCANSS. This implies that each combination of values
% determied by VARIATE, will be sampled SCANREP times.
% Example:
% SCANREP 100
% 
% NOTSIM := Prevents from doing any new stochastic simulation. It will
% make the corresponding plots though, by loading already existing
% outputs. Useful if you wish to change observables or add scatter
% plots without rerunning simulations. When used that way, time points
% should not be modified! Optional.
% Example: NOTSIM 1
%
% [] = parsejob(FILE)
%
% For more details see the documentation document in doc/ folder
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    ext = '\.job$'; % Valid extention for job files
    vsep = '\s+';

    dodot = 0;
    dodet = 0;
    doscandet = 0;
    doss = 0;
    dost = 0;
    doscanss = 0;
    notsim = 0;
    outfolder = '.';
    observe = {};
    scatter = cell(100,1);
    nscatter = 0;
    variate = cell(100,1);
    varvals = cell(100,1);
    nvals = zeros(100,1);
    nvariate = 0;
    name = '';

    if (~ regmbool(file, ext)) % Check extension
        status = 0;
    else
        [fid,mssg] = fopen(file); % Open file
        if (mssg) % Could not open the file...
            status = 0;
        else
            cprintf([1,0.5,0],'------');
            cprintf('text',' Parsing file %s ',file);
            cprintf([1,0.5,0],'------\n');
            lline = fgetl(fid);
            while ischar(lline) % For each line
                lline = regexprep(lline,'#.*','');
                vals = mysplit(lline,vsep);
                if (length(vals)>1)
                    sym = vals{1};
                    val = vals{2};
                    if(strcmp(sym,'') || strcmp(val,''))
                        % empty symbol or value
                    elseif(regmbool(sym,'^(#|%)'))
                        % comments
                    elseif(strcmp(sym,'NAME'))
                        name = val;
                    elseif(strcmp(sym,'MODEL'))
                        if(regmbool(val,'\.rxnm$'))
                            model = parserxnm(val);
                        elseif(regmbool(val,'\.sbproj$'))
                            sbproj = sbioloadproject(val);
                            model = simbio2model(sbproj.m1);
                        else
                            error('Invalid model extension');
                        end
                    elseif(strcmp(sym,'OUTFOLDER'))
                        outfolder = val;
                    elseif(strcmp(sym,'DOT'))
                        if(str2num(val))
                            dodot = 1;
                        end
                    elseif(strcmp(sym,'NOTSIM'))
                        if(str2num(val))
                            notsim = 1;
                        end
                    elseif(strcmp(sym,'SETINIT'))
                        setp = mysplit(val,':');
                        model = setinit(model,setp{1},str2num(setp{2}));
                    elseif(strcmp(sym,'SETDIST'))
                        setp = mysplit(val,':');
                        if(regmbool(setp{1},'^\/.*\/$'))
                            syms = cellregexm(model.ind2spc,setp{1}(2:end-1));
                        else
                            syms = mysplit(setp{1},',');
                        end
                        dist = load_table(setp{2});
                        model = addrvinit(model,syms,dist);
                    elseif(strcmp(sym,'SETCORR'))
                        setp = mysplit(val,':');
                        if(regmbool(setp{1},'^\/.*\/$'))
                            syms = cellregexm(model.ind2spc,setp{1}(2:end-1));
                        else
                            syms = mysplit(setp{1},',');
                        end
                        c = str2num(setp{2});
                        model = addcorr(model,syms,c);
                    elseif(strcmp(sym,'TYPE'))
                        types = mysplit(val,',|&');
                        for i = 1:length(types)
                            type = types{i};
                            if(strcmp(type,'DET'))
                                dodet = 1;
                            elseif(strcmp(type,'SS'))
                                doss = 1;
                            elseif(strcmp(type,'ST'))
                                dost = 1;
                            elseif(strcmp(type,'SCANDET'))
                                doscandet = 1;
                            elseif(strcmp(type,'SCANSS'))
                                doscanss = 1;
                            else
                                error('Invalid job type!');
                            end
                        end
                    elseif(strcmp(sym,'OBSERVE'))
                        if(strcmp(val,'ALL'))
                            observe = model.ind2spc;
                        else
                            observe = mysplit(val,',|&');
                        end
                    elseif(strcmp(sym,'DETPOINTS'))
                        detdt = eval(['[',val,']']);
                    elseif(strcmp(sym,'SSNREP'))
                        ssnrep = str2num(val);
                    elseif(strcmp(sym,'SSTPOINTS'))
                        ssts = eval(['[',val,']']);
                    elseif(strcmp(sym,'SCATTER'))
                        nscatter = nscatter + 1;
                        vals = mysplit(val,',|&');
                        scatter{nscatter} = vals;
                    elseif(strcmp(sym,'STNREP'))
                        stnrep = str2num(val);
                    elseif(strcmp(sym,'STTPOINTS'))
                        stts = eval(['[',val,']']);
                    elseif(strcmp(sym,'VARIATE'))
                        nvariate = nvariate + 1;
                        var = val;
                        varval = eval(['[',vals{3},']']);
                        if(length(varval)==1)
                            error(['VARIATE ',var,' has only one value. Use SETINIT']);
                        end
                        variate{nvariate} = var;
                        varvals{nvariate} = varval;
                        nvals(nvariate) = length(varvals{nvariate});
                    elseif(strcmp(sym,'SCANTIME'))
                        scantime = str2num(val);
                    elseif(strcmp(sym,'SCANREP'))
                        scanrep = str2num(val);
                    else
                        error(['Invalid job attribute option: ',sym]);
                    end
                end
                lline = fgetl(fid);
            end
        end
        
        fclose(fid);
        currentdir = cd; % save current folder
        s = 0;
        if(~(exist(outfolder,'dir'))) % we need to make the folder
                                      % first
            [s,r] = system(['mkdir ',outfolder]);
        end
        if(~(s==0))
            error('Outfolder does not exists and could not be created');
        else
            if(~regmbool(outfolder,'\/$'))
                outfolder = [outfolder,'/'];
            end
            cd(outfolder); % move to that folder
            cprintf([1,0.5,0],'------');
            cprintf('text',' Moving to %s ',outfolder);
            cprintf([1,0.5,0],'------\n');
            addpath(cd); % add path so that the function is
                         % recognized
            
            %% Save the model in rxnm format
            cprintf([1,0.5,0],'------');
            cprintf('text',' Saving Model ');
            cprintf([1,0.5,0],'------\n');
            model2rxnm(model);
            
            %% Make dot files
            if(dodot)
                cprintf([1,0.5,0],'------');
                cprintf('text',' Making Dot Graph ');
                cprintf([1,0.5,0],'------\n');
                dotfile = model2dot(model);
                png = dot2png(dotfile);
                img = imread(png);
                figure();
                image(img);
                axis off; % Remove axis ticks and numbers
                axis image; % Set aspect ratio to obtain square pixels 
            end
            
            % remove stuff that can not be used
            model = rmnonused(model);
            
            fname = [model.name,'_',name];
            
            %% Make deterministic calculations
            if(dodet)
                % TODO! do the same rate checking we do for stoch but for
                % deterministic
                if(~exist('detdt','var'))
                    error('DETPOINTS non specified');
                end
                cprintf([1,0.5,0],'------');
                cprintf('text',' Making Deterministic Simulation ');
                cprintf([1,0.5,0],'------\n');
                if(~(detdt(1)==0))
                    detdt = [0,detdt];
                end
                [Tdet,Xdet] = rundet(model,detdt);
                xtitles = model.ind2spc;
                for i = 1:length(xtitles)
                    xtitles(i) = regexprep(xtitles(i),'_','\\_');
                end
                XIDs = model.ind2spc;
                PIDs = model.ind2par;
                [X0,Params] = getinitvec(model);
                save([fname,'_det.mat'],'-mat','Tdet','Xdet','XIDs','X0','Params','PIDs');
                if(length(observe)>0)
                    obs = getobs(model,observe);
                    fig = figure(); hold on;
                    %set(fig,'Visible','Off');
                    nobs = length(obs);
                    colors = hsv(nobs);
                    for i = 1:nobs
                        subplot(nobs,1,i); hold on;
                        j = obs(i);
                        plot(Tdet,Xdet(:,j),'LineWidth',2,'Color',colors(i,:));
                        set(gca,'FontSize',12);
                        if(max(Xdet(:,j))>0)
                            axis([0 Tdet(end) 0 max(Xdet(:,j))]);
                        else
                            axis([0 Tdet(end) 0 1]);
                        end
                        xlabel('time','FontSize',20); ylabel(xtitles{j},'FontSize',20);
                    end
                    saveas(fig,[fname,'_det'], 'fig');
                    % add here to save figure
                end
            end
            
            if(doss || dost)
                % For stochastic simulations we need to check that rates makes sense
                % in the discrete and propenseties world
                [smodel] = cmstoch(model);
            end
            
            %% Time trajectories
            if(dost)
                if(~exist('stts','var') || ~exist('stnrep','var'))
                    error('STTPOINTS or STNREP non specified');
                end
                if(~(stts(1)==0))
                    stts = [0,stts];
                end
                simo = makesimo(smodel,name,'st',stts,stnrep);
                if(~notsim)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Making Moments Time Curves ');
                    cprintf([1,0.5,0],'------\n');
                    runsim(simo);
                    XIDs = smodel.ind2spc;
                    PIDs = model.ind2par;
                    [X0,Params] = getinitvec(model);
                    [Tst,Xst,STDst] = loadsimout(simo,XIDs);
                    % reshape data into a 3 dimensional matrix
                    X3st = zeros(length(Tst),length(XIDs),stnrep);
                    STD3st = zeros(length(Tst),length(XIDs),stnrep);
                    for i = 1:length(XIDs)
                        for n = 1:stnrep
                        X3st(:,i,n) = Xst{i}(:,n);
                        STD3st(:,i,n) = STDst{i}(:,n);
                        end
                    end
                    save([fname,'_st.mat'],'-mat','Tst','Xst', ...
                         'X3st','STDst','STD3st','XIDs','X0','Params','PIDs');
                end
                if(length(observe)>0)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Plotting observables trajectories ');
                    cprintf([1,0.5,0],'------\n');
                    fig = plotssim(simo,observe);
                end
            end
            
            %% States distributions
            if(doss)
                if(~exist('ssts','var') || ~exist('ssnrep','var'))
                    error('SSTPOINTS or SSNREP non specified');
                end
                if(~(ssts(1)==0))
                    ssts = [0,ssts];
                end
                simo = makesimo(smodel,name,'ss',ssts,ssnrep);
                if(~notsim)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Making States Distributions ');
                    cprintf([1,0.5,0],'------\n');
                    XIDs = smodel.ind2spc;
                    PIDs = model.ind2par;
                    [X0,Params] = getinitvec(model);
                    runsim(simo);
                    [Tss,Xss] = loadsimout(simo,XIDs);
                    % reshape data into a 3 dimensional matrix
                    X3ss = zeros(length(ssts),length(XIDs),ssnrep);
                    for i = 1:length(XIDs)
                        for n = 1:ssnrep
                            X3ss(:,i,n) = Xss{i}(:,n);
                        end
                    end
                    save([fname,'_ss.mat'],'-mat','Tss','Xss','X3ss','XIDs','X0','Params','PIDs');
                end
                if(length(observe)>0)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Plotting observables distribution ');
                    cprintf([1,0.5,0],'------\n');
                    fig = plotssim(simo,observe);
                end
                if(nscatter)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Making Scatter Plots ');
                    cprintf([1,0.5,0],'------\n');
                    scatter = scatter(1:nscatter);
                    plotscatter(simo,scatter);
                end
            end
            %% Scan values
            if(nvariate && (doscandet || doscanss))
                if(~exist('scantime','var'))
                    error('SCANTIME is not defined');
                end
                if(~(scantime(1)==0))
                    scantime = [0,scantime];
                end
                variate = variate(1:nvariate);
                varvals = varvals(1:nvariate);
                ntimes = length(scantime);
                nvals = nvals(1:nvariate);
                nscans = prod(nvals);
                [scans,indexes] = combsets(varvals);
                if(doscandet)
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Making Deterministic Scan ');
                    cprintf([1,0.5,0],'------\n');
                    XIDs = model.ind2spc;
                    PIDs = model.ind2par;
                    [DefX0,DefParams] = getinitvec(model);
                    XscanDET = zeros(nscans,length(XIDs),ntimes);
                    %tbins = 1:(scantime/1000):scantime;
                    smodel = model;
                    for s = 1:nscans
                        for v = 1:nvariate
                            smodel = setinit(smodel,variate{v},scans(s,v));
                        end
                        [tp,x] = rundet(smodel,scantime);
                        for i = 1:ntimes
                            ind = find(tp==scantime(i));
                            XscanDET(s,:,i) = x(ind,:)';
                        end
                    end
                    save([fname,'_scandet.mat'],'-mat','scantime','scans','indexes','variate',...
                         'varvals','XscanDET','XIDs','DefX0','DefParams','PIDs');
                    % plot only if we are variating 2 or 3 things
                    if((length(observe)>0) && (nvariate<3))
                        obs = getobs(model,observe);
                        nobs = length(obs);
                        colors = hsv(nobs);
                        for t = 1:(ntimes)
                            fig = figure(); hold on;
                            %set(fig,'Visible','Off');
                            for i = 1:nobs
                                subplot(nobs,1,i); hold on;
                                j = obs(i);
                                if(nvariate==1)
                                    plot(scans,XscanDET(:,j,t),'+--','LineWidth',2,'Color',colors(i,:));
                                    set(gca,'FontSize',12);
                                    axis([min(scans) max(scans) 0 (max(XscanDET(:,j,t))+1)]);
                                    xlabel(variate{1},'FontSize',20); ylabel(XIDs{j},'FontSize',20);
                                end
                                if(nvariate==2)
                                    Z = zeros(nvals(1),nvals(2));
                                    for s = 1:nscans
                                        x1 = indexes(s,1);
                                        x2 = indexes(s,2);
                                        Z(x1,x2) = XscanDET(s,j,t);
                                    end
                                    surf(varvals{1},varvals{2},Z');
                                    %scatter3(scans(:,1),scans(:,2),XscanDET(:,j));
                                    view(-60,60);
                                    set(gca,'FontSize',15);
                                    xlabel(variate{1},'FontSize',20); ylabel(variate{2},'FontSize',20);
                                    zlabel(XIDs{j},'FontSize',20);
                                    axis([min(scans(:,1)) max(scans(:,1)) min(scans(:,2)) max(scans(:,2)) 0 (max(XscanDET(:,j,t))+1)]);
                                end
                            end
                            set(fig,'name',['t=',num2str(scantime(t))]);
                            saveas(fig,[fname,'_scandet_',num2str(t)], 'fig');
                            % add here to save figure
                        end
                    end
                end
                if(doscanss)
                    if(~exist('scanrep','var'))
                        error('SCANREP is not defined');
                    end
                    cprintf([1,0.5,0],'------');
                    cprintf('text',' Making Stochastic Scan');
                    cprintf([1,0.5,0],'------\n');
                    XIDs = model.ind2spc;
                    PIDs = model.ind2par;
                    [~,DefParams] = getinitvec(model);
                    [smodel] = cmstoch(model);
                    XscanSS = cell(nscans,length(XIDs),ntimes);
                    X4scanSS = zeros(nscans,length(XIDs),ntimes,scanrep);
                    X3scanSS = zeros(nscans*scanrep,length(XIDs),ntimes);
                    scans3 = zeros(nscans*scanrep,nvariate);
                    indexes3 = zeros(nscans*scanrep,nvariate);
                    for s = 1:nscans
                        for v = 1:nvariate
                            smodel = setinit(smodel,variate{v},scans(s,v));  
                        end
                        simo = makesimo(smodel,[name,'_scan'],'ss',scantime,scanrep);
                        runsim(simo);
                        [~,ssX] = loadsimout(simo,XIDs);
                        for j = 1:length(XIDs)
                            X4scanSS(s,j,:,:) = ssX{j};
                            for i = 1:ntimes
                                XscanSS{s,j,i} = ssX{j}(i,:);
                            end
                        end
                    end
                    for s = 1:nscans
                        for n = 1:scanrep
                            k = (s-1)*scanrep + n;
                            X3scanSS(k,:,:) = X4scanSS(s,:,:,n);
                            scans3(k,:) = scans(s,:);
                            indexes3(k,:) = indexes(s,:);
                        end
                    end
                    save([fname,'_scanss.mat'],'-mat','scantime', ...
                         'scans','scans3','indexes','indexes3','scans3', ...
                         'variate','varvals', 'XscanSS', ...
                         'X3scanSS','X4scanSS','XIDs', 'DefParams','PIDs');
                    if((length(observe)>0) && (nvariate<3))
                        obs = getobs(smodel,observe);
                        nobs = length(obs);
                        colors = hsv(nobs);
                        for t = 1:ntimes
                            fig = figure(); hold on;
                            %set(fig,'Visible','Off');
                            for i = 1:nobs
                                j = obs(i);
                                means = zeros(nscans,1);
                                stds = zeros(nscans,1);
                                for s = 1:nscans
                                    means(s) = mean(XscanSS{s,j,t});
                                    stds(s) = std(XscanSS{s,j,t});
                                end
                                cvs = stds./means;
                                if(nvariate==1)
                                    subplot(nobs,2,(2*i-1)); hold on;
                                    plot(scans,means,'+--','LineWidth',2,'Color',colors(i,:));
                                    set(gca,'FontSize',12);
                                    axis([min(scans) max(scans) 0 (max(means)+1)]);
                                    xlabel(variate{1},'FontSize',20); ylabel(['<',XIDs{j},'>'],'FontSize',20);
                                    subplot(nobs,2,2*i); hold on;
                                    plot(scans,cvs,'+--','LineWidth',2,'Color',colors(i,:));
                                    set(gca,'FontSize',12);
                                    axis([min(scans) max(scans) 0 (max(cvs)+0.05*max(cvs))]);
                                    xlabel(variate{1},'FontSize',20); ylabel(['CV ',XIDs{j}],'FontSize',20);
                                elseif(nvariate==2)
                                    subplot(nobs,2,(2*i-1)); hold on;
                                    Z = zeros(nvals(1),nvals(2));
                                    for s = 1:nscans
                                        x1 = indexes(s,1);
                                        x2 = indexes(s,2);
                                        Z(x1,x2) = means(s);
                                    end
                                    surf(varvals{1},varvals{2},Z');
                                    %scatter3(scans(:,1),scans(:,2),means);
                                    view(-60,60);
                                    set(gca,'FontSize',12);
                                    xlabel(variate{1},'FontSize',20); ylabel(variate{2},'FontSize',20);
                                    zlabel(['<',XIDs{j},'>'],'FontSize',20);
                                    axis([min(scans(:,1)) max(scans(:,1)) min(scans(:,2)) max(scans(:,2)) 0 (max(means)+1)]);
                                    subplot(nobs,2,2*i); hold on;
                                    Z = zeros(nvals(1),nvals(2));
                                    for s = 1:nscans
                                        x1 = indexes(s,1);
                                        x2 = indexes(s,2);
                                        Z(x1,x2) = cvs(s);
                                    end
                                    surf(varvals{1},varvals{2},Z');
                                    %scatter3(scans(:,1),scans(:,2),cvs);
                                    view(-60,60);
                                    set(gca,'FontSize',12);
                                    maxcv = max([(max(cvs)+0.05*max(cvs)) 1]);
                                    axis([min(scans(:,1)) max(scans(:,1)) min(scans(:,2)) max(scans(:,2)) 0 maxcv]);
                                    xlabel(variate{1},'FontSize',20); ylabel(variate{2},'FontSize',20);
                                    zlabel(['CV ',XIDs{j}],'FontSize',20);
                                end
                            end
                            set(fig,'name',['t=',num2str(scantime(t))]);
                            saveas(fig,[fname,'_scanss_',num2str(t)], 'fig');
                        end
                        if(nscatter)
                            %cprintf([1,0.5,0],'------');
                            %cprintf('text',' Making Scatter Scans Plots ');
                            %cprintf([1,0.5,0],'------\n');
                        end
                    end
                end
            end
            cd(currentdir) % move back to the previous folder
        end
    end
end