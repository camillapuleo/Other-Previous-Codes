close all
clear
clear classes % update classes in the memory
addpath('utils');



%% Flat simulated consumption path using migration model
m4=my_model_migration;
m4.ngridm=500;
m4.df=1/(1+m4.r); %flat consumption hopefully
m4.sigma=0;
m4.lambda=eps; %no EV taste shocks
m4.nsims=5;
m4.init=[5 20];
tic
m4.solve_dcegm;
t=toc;
fprintf('Migration model solved with DC-EGM in %s\n',ht(t));
m4.plot('policy')
m4.plot('value')
m4.plot('prob_stay')
m4.sim;
m4.plot('sim consumption');

fprintf('\nPress any key to continue..');pause;fprintf('\n\n');

%% Nice simulation graphics using migration model
m5=my_model_migration;
m5.ngridm=500;
m5.df=1/(1+m5.r); %flat consumption hopefully
m5.sigma=0.35;
m5.lambda=0.2; %some EV taste shocks
m5.nsims=50;
m5.init=[5 20];
tic
m5.solve_dcegm;
t=toc;
fprintf('Migration model solved with DC-EGM in %s\n',ht(t));
m5.plot('policy');
m5.plot('value');
m5.plot('prob_stay');
m5.sim;
m5.plot('sim');
fprintf('Simulation plots for migration model produced\n')

m5.plot('prob_work');
m5.sim;
m5.plot('sim');
fprintf('Simulation plots for retirement model produced\n')
