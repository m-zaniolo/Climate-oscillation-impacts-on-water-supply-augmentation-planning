function sim_param = generate_sim_param(inflow, osc_short,osc_med,osc_long, storage, N, M, k, slope, intercept)

sim_param.inflow        = inflow/12; %MCM/month
sim_param.storage       = storage;
sim_param.min_release   = 0;
sim_param.lsv           = 0;
sim_param.s0            = 15;
sim_param.H             = 1200; % length of simulation horizon in months
sim_param.Ny            = 100;
sim_param.T             = 12;
sim_param.nsim          = size(inflow, 1);
sim_param.hedge_param   = 15;
sim_param.water_demand  = 90/12; %MCM/month
sim_param.osc_short     = osc_short;
sim_param.osc_med       = osc_med;
sim_param.osc_long      = osc_long;
%sim_param.operation_param = [0.4, 0.4, 0.4, 0.4]; %temporary

%% tech
%damcost_model              = load('data/cost_model_exp.mat');
%idcost                     = find(damcost_model.cost_model_exp.storage == storage);
%sim_param.damcost          = damcost_model.cost_model_exp.dam_cost(idcost);
%sim_param.storage_exp_cost = damcost_model.costmodel.dam_cost(1); % expanding storage means creating a new 50 MCM storage
%% tech parameters 
sim_param.tech_cap    = 1;

sim_param.intercept = intercept;
sim_param.slope = slope;

%% optimization parameters
sim_param.N                = N; 
sim_param.M                = M; 
sim_param.k                = k; 
