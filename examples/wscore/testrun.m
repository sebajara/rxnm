clear all; close all;

addpaths();

model = parserxnm('wscore.rxnm');
ts = [0,100000];
s0_vec = logspace(0,5,10);

%% Normal MATLAB run
model2der(model); % Make derivative function
model2eval(model); % Make eval function
tmodel = model;
lam1 = [];
tic;
for i = 1:length(s0_vec)
    tmodel = setinit(tmodel,'s0',s0_vec(i));
    [~,x]  = rundet(tmodel,ts);
    lam1(i) = evalinmodel(tmodel,ts(end),x(end,:),'lam');
end
toc;

%% ODEMEX run
model2der(model,'Mode',1,'Compile',1);
%model2der(model,'Mode',1);
tmodel = model;
lam2 = [];
tic;
for i = 1:length(s0_vec)
    tmodel = setinit(tmodel,'s0',s0_vec(i));
    [~,x]  = rundet(tmodel,ts,'Mode',1);
    lam2(i) = evalinmodel(tmodel,ts(end),x(end,:),'lam');
end
toc;

figure();
subplot(1,2,1);
plot(log(s0_vec),lam1,'--o');
subplot(1,2,2);
plot(log(s0_vec),lam2,'--o');
