function [T,X,STD] = loadsimoout (simo,obscel)

% This function takes a SIM structure and a string cell array
% containing species names (OBSCEL). Assuming that the simulation was
% 'executed', it will parse the output files and obtain the
% corresponding data for the species specified in OBSCEL.
% 
% For 'ss' simulations, T will be a 1,1 cell array, X a NOBS,1 cell array,
% each containing a NTIMEPOINTS by NREP matrix with the state of
% the specie at that time. STD is not instantiated.
% 
% For 'st' simulations, T will be a 2,1 cell array, containing the upper and
% lower values of the time intervals. X is a NOBS,1 cell array, each
% containing a NTIMEPOINTS by NREP matrix with the average value of
% the specie during that time interval. STD is the same size, but
% containing the standard deviation during that interval.
% 
% [T,X,STD] = loadsimout(SIMO,OBSCEL)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
name = simo.name;
model = simo.model;
type = simo.type;
ts = simo.ts;
nrep = simo.nrep;

% Get the index of each observable
observe = getobs(model,obscel);
nobs = length(observe);

species = fieldnames(model.species);
nspec = length(species);
ntps = length(ts);

if (strcmp(type,'ss'))
    % convention for the output
    filestr = [model.name,'_',name,'_ss'];
    seriestr = '%f';
    for i = 1:nspec
        seriestr = [seriestr, ' %f'];
    end
    file = [filestr,'.out'];
    fid = fopen(char(file));
    data = textscan(fid,seriestr);
    fclose(fid);
    T = data{1}(1:ntps);
    X = cell(nobs,1);
    for obs = 1:nobs
        X{obs,1} = zeros(ntps,nrep);
        o = observe(obs);
        j = o + 1;
        for n = 1:nrep
            tps = [(ntps*(n-1)+1):ntps*n];
            x = data{1,j}(tps);
            X{obs,1}(:,n) = x;
        end
    end
elseif(strcmp(type,'st'))
    % convention of for the output
    filestr = [model.name,'_',name,'_st'];
    seriestr = '%f %f';
    for i = 1:(2*nspec)
        seriestr = [seriestr, ' %f'];
    end
    X = cell(nobs,1);
    STD = cell(nobs,1);
    for obs = 1:nobs
        X{obs,1} = zeros(ntps,nrep);
        STD{obs,1} = zeros(ntps,nrep);
    end
    for r = 1:nrep
        file = [filestr,num2str(r-1),'.out'];
        fid = fopen(char(file));
        data = textscan(fid,seriestr);
        fclose(fid);
        T = zeros(ntps,2);
        if(length(data{1})<ntps)
            error('Less time points than expected');
        end
        T(:,1) = data{1};
        T(:,2) = data{2};
        for obs = 1:nobs
            o = observe(obs);
            j = o + 2;
            x = data{j};
            x2 = data{j+nspec};
            std = sqrt(x2 - x.^2);
            X{obs,1}(:,r) = x;
            STD{obs,1}(:,r) = std;
        end
    end
end

clearvars data;

end