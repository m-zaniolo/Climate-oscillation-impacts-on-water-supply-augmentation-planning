function J = obj_func_planning( param, sim_param )

inflow               = sim_param.inflow;
s0                   = sim_param.s0;
H                    = sim_param.H;
Ny                   = sim_param.Ny;
T                    = sim_param.T;
nsim                 = sim_param.nsim;
osc_short            = sim_param.osc_short;
osc_med              = sim_param.osc_med;
osc_long             = sim_param.osc_long;
w                    = sim_param.water_demand;
hedge_param          = sim_param.hedge_param;
N                    = sim_param.N;
M                    = sim_param.M;
K                    = sim_param.k;
tech_cap             = sim_param.tech_cap;
storage              = sim_param.storage;
intercept            = sim_param.intercept;
slope                = sim_param.slope;

%% simulation
irr_cost_tot = nan(1, nsim);
plan_costs_tot = nan(1, nsim);

for sim = 1:nsim
    % pre-allocate vectors
    q = [0, inflow(sim, :)];
    short = osc_short(sim, :);
    med = osc_med(sim, :);
    long = osc_long(sim, :);
    % initial conditions
    s_ = s0;
    irr_cost = 0;
    planning_costs  = 0; 
    capacity    = zeros(1,H);
    
    for t = 1:H
        w_ = w-capacity(t); 

        %reservoir mass balance
        hedge = hedge_param * s_;
        flood = 10*s_*(1/2 - w_/storage) - 9*(storage/2 - w_) + w_ ;
        u = max( min(hedge, w_), flood);
        [s_, r] = mass_balance(s_, u, q(t+1), storage);
        irr_def = max(0, w_-r);

        if t>120 %spin up time
            irr_cost = irr_cost + irr_def^2;
        end
        
        % evaluate input values
        norm_inputs = [s_/storage, capacity(t)/6];
        if t>6
            if ~isnan(short(t))
                delta = short(t) - short(t-6) ;
                norm_inputs = [norm_inputs, short(t), delta];
            end
            if ~isnan(med(t))
                delta = med(t) - med(t-6) ;
                norm_inputs = [norm_inputs, med(t), delta];
            end
            if ~isnan(long(t))
                delta = long(t) - long(t-6) ;
                norm_inputs = [norm_inputs, long(t), delta];
            end
        else
            if ~isnan(short(t))
                norm_inputs = [norm_inputs, short(t), 0];
            end
            if ~isnan(med(t))
                norm_inputs = [norm_inputs, med(t), 0];
            end
            if ~isnan(long(t))
                norm_inputs = [norm_inputs, long(t), 0];
            end

        end

        % evaluate planning policy
        y = rbf(param, N , M, K, norm_inputs);
        tf = min(t+13, H);
        capacity(t+1:tf) = floor(y*7); %gives you capacity to be installed 
        

    end
    irr_cost_tot(sim)     = irr_cost/90; %divided by number of years in sim horizon minus spin up
    costs                 = capacity2costs(capacity(120:end), slope, intercept) ;
    plan_costs_tot(sim)   = sum(costs);

end

infr_cost = plan_costs_tot; 
J = [mean(irr_cost_tot), mean(infr_cost)];

end
