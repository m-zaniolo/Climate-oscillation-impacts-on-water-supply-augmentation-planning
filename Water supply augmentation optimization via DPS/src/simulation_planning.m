function [Jirr, Jinf, s, r, irr_def, planning_costs, lifetimes, capacity, delta] = simulation_planning( param, sim_param, inflow, idx )

inflow               = inflow/12; % inflow is measured in MCM per year
s0                   = sim_param.s0;
H                    = sim_param.H;
osc_short            = sim_param.osc_short;
osc_med              = sim_param.osc_med;
osc_long             = sim_param.osc_long;
w                    = sim_param.water_demand;
hedge_param          = sim_param.hedge_param;
N                    = sim_param.N;
M                    = sim_param.M;
K                    = sim_param.k;
storage              = sim_param.storage;
intercept            = sim_param.intercept;
slope                = sim_param.slope;


%% simulation

s               = nan(size(inflow));   % storage
r               = nan(size(inflow));   % release
irr_def         = nan(size(inflow));   % irrigation deficit
irr_cost        = nan(size(inflow));   % irrigation cost
w_              = nan(size(inflow));   % residual water demand after technology production
capacity        = zeros(size(inflow)); % installed capacity


% pre-allocate vectors
q = [0, inflow];
short = osc_short(idx, :);
med = osc_med(idx, :);
long = osc_long(idx, :);

% initial conditions
s(1)        = s0;


for t = 1:H
    w_(t) = w - capacity(t);

    hedge = hedge_param * s(t);
    flood = 10*s(t)*(1/2 - w_(t)/storage) - 9*(storage/2 - w_(t)) + w_(t) ;
    u = max( min(hedge, w_(t)), flood);
    [s(t+1), r(t+1)] = mass_balance(s(t), u, q(t+1), storage);
    irr_def(t) = max(0, w_(t)-r(t+1));
    
    irr_cost(t) = irr_def(t)^2;

    % evaluate input values
    norm_inputs = [s(t+1)/storage, capacity(t)/6];
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
        capacity(t+1:tf) = floor(y*7); 
        
end
    % extract planning costs
    [costs, lifetimes] = capacity2costs(capacity(120:end), slope, intercept);
    planning_costs = sum(costs);

% calculate objectives
Jinf = mean(planning_costs) ;
irr_cost = irr_cost(:, 121:end); % irrigation costs calculated after spin up time
Jirr = sum(irr_cost)/90;


end
