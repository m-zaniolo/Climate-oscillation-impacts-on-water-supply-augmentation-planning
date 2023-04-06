clear
clc
addpath('src/')
addpath('configuration_files/')


%% run for baseline case
conf = 'conf_l_d30_s0.5_i2';

%% read configuration file
[ osci_type, dam_size, slope, intercept] = import_conf( strcat(conf, '.txt') );

%% set optimization param
% parameters are set low for demo optimization. Results won't be converged.
nsim        = 1;   % number of climate scenarios
Niter       = 1;   % number of seeds for GA
PopSize     = 30;  % GA population size
numGen      = 100; % GA number of generations

% optim param used in large scale optimization for paper results:
%nsim        = 10;  % number of climate scenarios
%Niter       = 20;  % number of seeds for GA
%PopSize     = 80;  % GA population size
%numGen      = 300; % GA number of generations


k = 1; %output
M = 2; %inputs: storage, capacity



%% load variables for objectives computation
global sim_param

% load inflow
inflow    = load( strcat('data/runoff_', osci_type, '.csv'));% runoff_', inflow_type, '.csv') );

% load oscillation signals
if exist(strcat('data/osi_', osci_type, '_s.csv'), "file")>0
    osc_short = load( strcat('data/osi_', osci_type, '_s.csv'))'; 
    M = M+2; % include value and delta of oscillation signal in the input set
    osc_short = osc_short/2+0.5;
    m_osc_short = nan(size(osc_short,1), 1200);
    for i=1:size(osc_short,1)
        m_osc_short(i,:) = [repmat( osc_short(i,1), 1, 6) ,  [interp1qr(6:12:1200,osc_short(i,:),  7:1200-6)]', repmat( osc_short(i,end), 1, 6)];
    end
else
    m_osc_short = NaN(100,1200);
end

if exist(strcat('data/osi_', osci_type, '_m.csv'), "file")>0
    osc_med = load( strcat('data/osi_', osci_type, '_m.csv'))'; 
    M = M+2; % include value and delta of oscillation signal in the input set
    osc_med   = osc_med/2+0.5;
    m_osc_med = nan(size(osc_med,1), 1200);
    for i=1:size(osc_med,1)
        m_osc_med(i,:) = [repmat( osc_med(i,1), 1, 6) ,  [interp1qr(6:12:1200,osc_med(i,:),  7:1200-6)]', repmat( osc_med(i,end), 1, 6)];
    end

else
    m_osc_med = NaN(100,1200);
end

if exist(strcat('data/osi_', osci_type, '_l.csv'))>0
    osc_long = load( strcat('data/osi_', osci_type, '_l.csv'))'; 
    M = M+2; % include value and delta of oscillation signal in the input set
    osc_long  = osc_long/2+0.5;
    m_osc_long = nan(size(osc_long,1), 1200);
    for i=1:size(osc_long,1)
        m_osc_long(i,:) = [repmat( osc_long(i,1), 1, 6) ,  [interp1qr(6:12:1200,osc_long(i,:),  7:1200-6)]', repmat( osc_long(i,end), 1, 6)];
    end
else
    m_osc_long = NaN(100,1200);
end

N = M+k+1; % number of basis is set following rule of thumb: #inputs + #outputs + 1

%% set simulation parameters

sim_param  = generate_sim_param(inflow(1:nsim, :),m_osc_short(1:nsim, :), m_osc_med(1:nsim, :), m_osc_long(1:nsim, :),...
    dam_size, N, M, k, slope, intercept);

%% set optimization parameters

n_param = N*(2*M + k) + k; 
lb_node = [repmat([-1, 0], 1, M), zeros(1,k)] ; % lower bounds of node param
lb = [zeros(1,k), repmat(lb_node, 1, N) ]; % lower bound of RBF param
ub = ones(1,n_param); % upper bound of RBF param
kk = []; % initialize parameters 
JJ = []; % initialize objectives


%%
for iter = 1:Niter
    FitnessFunction = @(x) obj_func_planning(x, sim_param);
    options = optimoptions('gamultiobj', 'UseParallel', false, 'PopulationSize', PopSize, 'Generations', numGen); 
    [k, J] = gamultiobj(FitnessFunction, n_param,[],[],[],[],lb,ub, options); 
    kk = [kk;k];
    JJ = [JJ;J];
    
end

%% identify set of pareto efficient solutions 

[~, idx_un] = unique(JJ,'rows');

Ju = JJ(idx_un, :);
ku = kk(idx_un, :);

indPar=non_dominated_front(Ju');
J = Ju(indPar,:);
k = ku(indPar,:);

[~, index_filter] = filter_reference(J, [3, 3]);
Jfilter = J(index_filter,:);
kfilter = k(index_filter,:);

%% simulate optimal solutions
for i = 1:size(kfilter,1)
    [Jirr(i), Jinf(i), s(i,:), r(i,:), irr_def(i,:), planning_costs(i,:), avg_lifetime{i}, capacity(i,:)]  = simulation_planning( kfilter(i,:), sim_param, inflow(1, :) , 1);
end


%% plots
figure; scatter(J(:,1), J(:,2), 'filled')
xlabel('annual squared irrigation deficit'); ylabel('infrastructure cost')
set(gca, 'FontSize', 12); box on; grid on;
title('Unconverged results for demo purpose')
srt_fig = strcat( 'figures/pf_', conf,'.png');
saveas(gcf,srt_fig)


sim = i; % displaying simulation for least deficit solution
figure;
t = datetime(2000, 1, 1)+calmonths(0:1200-1);
subplot(411); plot(t, m_osc_short(1, :), 'LineWidth', 1.5);
hold on; plot(t, m_osc_med(1, :), 'LineWidth', 1.5); 
plot(t, m_osc_long(1, :), 'LineWidth', 1.5);

str_title = strcat('Unconverged results for demo only' );
title(str_title)
ylabel('normalized oscillations'); set(gca, 'FontSize', 12);

subplot(413); plot(t, s(sim, 1:end-1), 'LineWidth', 1.5);
hold on; plot(t(1:1200), repmat(dam_size, 1200, 1), 'k--'); set(gca, 'FontSize', 12);
ylim([0, 60])
grid on; title('Storage'); ylabel('[MCM]')

subplot(412); %plot(t, r(sim, 1:end-1), 'LineWidth', 1.5);
infl = sum( reshape( inflow(sim, :)/12, 12, 100 ) );
plot(t(1:12:end), infl, 'LineWidth', 1.5);
hold on;  plot(t(1:1200), repmat(90, 1200, 1), 'k--');
hold on; plot(t, irr_def(sim, :), 'LineWidth', 1.5);

set(gca, 'FontSize', 12);
legend('inflow', 'demand', 'deficit')
grid on; title('Annual Inflow'); ylabel('[MCM]')


subplot(414); area(t, capacity(sim, 1:length(t)), 'LineWidth', 1.5);  grid on;
legend('capacity'); set(gca, 'FontSize', 12)
ylim([0,5]);
srt_fig = strcat( 'figures/traj_', conf,'.png');
saveas(gcf,srt_fig)

