%--------------------------------------------------------------------------
% Summarize the results of model fitting and Simulate the effects of
% interventions with the cooperative infection model
%--------------------------------------------------------------------------

% Add the folder 'Functions' to the MATLAB search path
addpath('Functions');

%% Reorganize figures for comparison of simulation and actual data
% Load data
Countries_left = {'France','Spain','Italy','Germany','UK','USA','Turkey'};
ActiveCases = readtable('ActiveCases.csv','ReadRowNames',true);
TotalCases = readtable('TotalCases.csv','ReadRowNames',true);
Countries = TotalCases.Properties.RowNames;
Dates = TotalCases.Properties.VariableNames;
Start_Date = 'Feb15';
date0_idx_common = find(strcmp(Dates,Start_Date));

% Define colormap and figure handles
colors = brewermap(8,'Set2');
raw_I_x = cell(1,7);
raw_I_y = cell(1,7);
sim_I_x = cell(1,7);
sim_I_y = cell(1,7);
sample_I_x = cell(1,7);
sample_I_y = cell(1,7);
raw_N_x = cell(1,7);
raw_N_y = cell(1,7);
sim_N_x = cell(1,7);
sim_N_y = cell(1,7);
sample_N_x = cell(1,7);
sample_N_y = cell(1,7);

% Read data from the original fig files
for i=1:length(Countries_left)
    fig_name=strcat(Countries_left{i},'-SIRD cooperative infection.fig');    
    idx=find(strcmp(Countries,Countries_left{i}));
    if ActiveCases{idx,1}>0
        start_idx = 1;
    else
        start_idx = max(find(ActiveCases{idx,:}==0))+1;
    end
    delta_t = start_idx - date0_idx_common;
    
    open(fig_name);
    h1=subplot(1,2,1);
    xd_handles = findobj(h1,'-property','xdata');
    raw_I_x{i} = xd_handles(1).XData+delta_t;
    raw_I_y{i} = xd_handles(1).YData;
    sim_I_x{i} = xd_handles(2).XData+delta_t;
    sim_I_y{i} = xd_handles(2).YData;
    sample_I_x{i} = xd_handles(3).XData+delta_t;
    sample_I_y{i} = xd_handles(3).YData;
    
    h2=subplot(1,2,2);
    xd_handles = findobj(h2,'-property','xdata');
    raw_N_x{i} = xd_handles(1).XData+delta_t;
    raw_N_y{i} = xd_handles(1).YData;
    sim_N_x{i} = xd_handles(2).XData+delta_t;
    sim_N_y{i} = xd_handles(2).YData;
    sample_N_x{i} = xd_handles(3).XData+delta_t;
    sample_N_y{i} = xd_handles(3).YData;
    
    close(gcf);
end

% Reorganize figures, merge the figures from different countries
figure_I = figure;
figure_N = figure;
for i=1:7
    figure(figure_I);
    subplot(2,4,i);
    hold on;box on;   
    fill(sample_I_x{i},sample_I_y{i},[0.8 0.8 0.8],'EdgeColor','none');    
    plot(sim_I_x{i},sim_I_y{i},'Color',colors(2,:));
    scatter(raw_I_x{i},raw_I_y{i},[],colors(1,:));
    xlabel('Days since Feb 15');
    ylabel('Number of people');
    xlim([0 80]);
    title(Countries_left{i});
    if i==7
        legend('Simulation (all sampled parameter sets)','Simulation (optimal fitting)',...
            'Actual data');
    end
    
    figure(figure_N);
    subplot(2,4,i);
    hold on;box on;    
    fill(sample_N_x{i},sample_N_y{i},[0.8 0.8 0.8],'EdgeColor','none');    
    plot(sim_N_x{i},sim_N_y{i},'Color',colors(2,:));
    scatter(raw_N_x{i},raw_N_y{i},[],colors(1,:));
    xlabel('Days since Feb 15, 2020');
    ylabel('Number of people');
    xlim([0 80]);
    title(Countries_left{i});
    if i==7
        legend('Simulation (all sampled parameter sets)','Simulation (optimal fitting)',...
            'Actual data');
    end
end

%% Load parameter estimation results for each country
Countries = {'USA','Italy','Spain','Germany','France','Iran',...
    'UK','Turkey','Switzerland'};
Models = {'SIRD original','SIRD virus pool','SIRD cooperative infection'};
ModelNames = {'Original SIRD','SIRD with virus pool','SIRD with virus pool and cooperative infection'};

Params_LBs = {[1e-6 1e-12],[1e-8 1 1e-8 1e-6],...
    [1e-8 1 1e-8 1e-6 1.1 1e-3]};
Params_UBs = {[1e2 1e-6],[1e-2 100 1e-2 1e2],...
    [1e-2 100 1e-2 1e2 5 1e2]};
ODE_Funcs = {@oriSIRD,@modSIRD_Vpool,@modSIRD_Vpool_CBind};
Jacobian_Funcs = {@oriSIRD_Jacobian,@modSIRD_Vpool_Jacobian,@modSIRD_Vpool_CBind_Jacobian};
Virus_Tags = [0 1 1];

params_samps_model3=[];
country_tags=[];
params_min_model3=[];
params_max_model3=[];
country_tags_unique=[];
for i_country=1:9
    for i_model=1:3
        run_name = strcat(Countries{i_country},'-',Models{i_model});
        load(strcat(run_name,'.mat'));
        Opt_Params{i_country,i_model}=10.^x_opt;
        Opt_Fobj(i_country,i_model)=f_opt/(size(ActiveCases,2)+1-start_idx);
        Sample_Params{i_country,i_model}=10.^x_kept;
        Sample_Fobj{i_country,i_model}=f_kept;
    end
    if size(x_kept,1)>1
        params_samps_model3=[params_samps_model3;Sample_Params{i_country,3}];
        country_tags=[country_tags;repmat(Countries(i_country),size(x_kept,1),1)];
        params_min_model3=[params_min_model3;min(Sample_Params{i_country,3})];
        params_max_model3=[params_max_model3;max(Sample_Params{i_country,3})];
        country_tags_unique=[country_tags_unique;Countries(i_country)];
    end
end

%% Plot optimal cost function values reached by each model
figure;
violinplot(log10(Opt_Fobj),ModelNames);xtickangle(45);ylim([-4 0]);
title('Distribution of optimal cost function values');
ylabel('log_{10}(Cost function)');
hold on;
plot([0.5 3.5],[-2.3010 -2.3010],':','Color',[0.5 0.5 0.5]);box on;

%% Comparing model parameters across countries
Model3_params_names = {'c (virus shedding)','d (virus decay)',...
    'k (susceptibility to infection)','a (patient clearance)','n (cooperativity in infection)',...
    'K (saturation in response to viral dose)'};
Model3_params_symbols = {'c','d','k','a','n','K'};

figure;
for i=1:6
    subplot(2,3,i);hold on;
    for j=1:length(country_tags_unique)
        fill([0.6 1.4 1.4 0.6]+(j-1)*[1 1 1 1],[log10(params_min_model3(j,i)) log10(params_min_model3(j,i)) ...
            log10(params_max_model3(j,i)) log10(params_max_model3(j,i))],[0.9 0.9 0.9]);
    end
    xticks(1:length(country_tags_unique));
    xticklabels(country_tags_unique);
    xtickangle(45);
    title(Model3_params_names{i});
    ylabel(strcat('log_{10}',Model3_params_symbols{i}));
    box on;
end

%% Simulate the cooperative infection term and effects of interventions
n_keep=5000; % Number of sampled parameter sets kept for the following analysis

% Simulate the cooperative infection term
Hill_res_1=zeros(n_keep,7); %Remaining hill term when virus pool is cut to 1% of the maximal size
Hill_res_10=zeros(n_keep,7); %Remaining hill term when virus pool is cut to 10% of the maximal size
colors = brewermap(8,'Set2');
figure;
Hill_res_1_10_opt=zeros(2,7);
for i = 1:length(country_tags_unique)
    i_country=find(strcmp(Countries,country_tags_unique{i}));
    country_tags_unique{i}
    
    if ActiveCases{i_country,1}>0
        start_idx = 1;
    else
        start_idx = max(find(ActiveCases{i_country,:}==0))+1;
    end
    
    data_cases_active = ActiveCases{i_country,start_idx:end};
    days_data_available = 0:length(data_cases_active)-1;
    
    start_date = Dates{start_idx};
    y_initial=[TotalPopulations(i_country)-TotalCases{i_country,start_idx} ActiveCases{i_country,start_idx} ...
        TotalCases{i_country,start_idx}-ActiveCases{i_country,start_idx} 0];
    
    % Retrieve the sampled parameter sets and optimal parameter set
    sample_params = Sample_Params{i_country,3};
    sample_fobj = Sample_Fobj{i_country,3};
    opt_params = Opt_Params{i_country,3};
    [~,idx_opt]=min(sample_fobj);
    nsample=size(sample_params,1);
    idx=randperm(nsample,n_keep);
    if ~ismember(idx_opt,idx) % Force the optimal solution to be kept
        idx(end)=idx_opt;
    end   
    sample_params=sample_params(idx,:);
    
    %Define the cooperative infection term
    V_mat=zeros(n_keep,101);
    Hill_mat=zeros(n_keep,101);
    
    % Create the subplot
    subplot(2,4,i);
    hold on;
    
    % Perform simulation for the optimal parameter set
    yo=SimDiseaseDynamics(opt_params,y_initial,size(ActiveCases,2)-start_idx);
    V_max=max(yo(:,4));
    n=opt_params(5);
    K=opt_params(6);
    V_opt=0:V_max/100:V_max;
    Hill_opt=(V_opt.^n)./(V_opt.^n+K);
    
    % Simulate with the sampled parameter sets
    for j=1:n_keep
        yo=SimDiseaseDynamics(sample_params(j,:),y_initial,size(ActiveCases,2)-start_idx);
        V_max=max(yo(:,4));
        n=sample_params(j,5);
        K=sample_params(j,6);
        V_mat(j,:)=0:V_max/100:V_max;
        Hill_mat(j,:)=(V_mat(j,:).^n)./(V_mat(j,:).^n+K);
    end
    
    % Compute cooperative infection term for 90% and 99% reduction of virus
    Hill_res_1(:,i)=Hill_mat(:,2)./max(Hill_mat,[],2);
    Hill_res_10(:,i)=Hill_mat(:,11)./max(Hill_mat,[],2);
    Hill_res_1_10_opt(1,i)=Hill_opt(2)/max(Hill_opt);
    Hill_res_1_10_opt(2,i)=Hill_opt(11)/max(Hill_opt);
    
    % Plot the results
    xx=0:0.01:1;
    fill([xx xx(end:-1:1)],[min(Hill_mat./max(Hill_mat,[],2)) max(Hill_mat(:,end:-1:1)./max(Hill_mat,[],2))],...
            [0.8 0.8 0.8],'EdgeColor',[1 1 1]);
    plot(xx,Hill_opt/max(Hill_opt),'Color',colors(1,:));
    plot(xx,median(Hill_mat./max(Hill_mat,[],2)),'Color',colors(2,:));
    title(country_tags_unique{i});
    xlabel('[V]/[V]_{max}');ylabel('H/H_{max}');
    legend('All sampled parameter sets','Parameter set with optimal fitting',...
        'Median values');
end
% Plot the results for 90% and 99% reduction of viral dose
figure;
subplot(1,2,1);
violinplot(Hill_res_10,country_tags_unique,'ShowData',false,'BandWidth',0.002);
hold on;
scatter(1:7,Hill_res_1_10_opt(2,:),150,colors(7,:),'p','filled','MarkerEdgeColor',[0 0 0]);
title('Viral dose reduced by 90%');
ylim([0 1]);xtickangle(45);ylabel('Relative cooperative infection term');
subplot(1,2,2);
violinplot(Hill_res_1,country_tags_unique,'ShowData',false,'BandWidth',0.002);
hold on;
scatter(1:7,Hill_res_1_10_opt(1,:),150,colors(7,:),'p','filled','MarkerEdgeColor',[0 0 0]);
title('Viral dose reduced by 99%');
xtickangle(45);ylabel('Relative cooperative infection term');

% Simulate different intervention strategies

% Define parameters
alpha=0.1;
day_max = 100;
day_start_intervention = 48;
t_all = 0:day_max;
for i=1:2
    f_I{i}=figure;
end

median_efficiency=zeros(7,6);
opt_efficiency=zeros(7,6);
for i = 1:length(country_tags_unique)
    i_country=find(strcmp(Countries,country_tags_unique{i}));
    country_tags_unique{i}
    
    if ActiveCases{i_country,1}>0
        start_idx = 1;
    else
        start_idx = find(ActiveCases{i_country,:}==0, 1, 'last' )+1;
    end
    delta_t = start_idx - 1;
    
    data_cases_active = ActiveCases{i_country,start_idx:end};
    days_data_available = 0:length(data_cases_active)-1;
    
    start_date = Dates{start_idx};
    y_initial=[TotalPopulations(i_country)-TotalCases{i_country,start_idx} ActiveCases{i_country,start_idx} ...
        TotalCases{i_country,start_idx}-ActiveCases{i_country,start_idx} 0];
    
    % Retrieve sampled parameter sets
    sample_params = Sample_Params{i_country,3};
    sample_fobj = Sample_Fobj{i_country,3};
    opt_params = Opt_Params{i_country,3};
    [~,idx_opt]=min(sample_fobj);
    nsample=size(sample_params,1);
    idx=randperm(nsample,n_keep);
    if ~ismember(idx_opt,idx) % Force the optimal solution to be kept
        idx(end)=idx_opt;
    end   
    sample_params=sample_params(idx,:);
        
    I_original = zeros(n_keep,day_max+1);
    I_intervention = zeros(n_keep,day_max+1);
    
    for j=1:n_keep
        yo=SimDiseaseDynamics(sample_params(j,:),y_initial,day_max);
        I_original(j,:)=yo(:,2)'; % Currently infected        
    end
    
    % Save the distributions of I and T on day_max
    I_end_sample=zeros(n_keep,4);    
    I_end_sample(:,1)=I_original(:,end);
        
    % Visualize results for original models (i.e. no intervention)
    for ni=1:2
        figure(f_I{ni});
        subplot(2,4,i);
        fill([t_all t_all(end:-1:1)]+delta_t,[min(I_original) max(I_original(:,end:-1:1))],...
            [0.8 0.8 0.8],'EdgeColor',[1 1 1]);
        hold on;
    end
    
    % Simulate and visualize interventions
    for ni = 1:2
        lists_alpha=[0.5 0.2 0.1];
        for k=1:3
            alpha=lists_alpha(k);
            for j=1:n_keep
                if ni == 1 % Population-based intervention
                    yi=SimIntervention(sample_params(j,:),...
                        [alpha 1 1 1 1 alpha^(-sample_params(j,5))],...
                        y_initial,day_start_intervention,day_max);
                else       % Infected individual-based intervention
                    yi=SimIntervention(sample_params(j,:),...
                        [alpha 1 1 1 1 1],...
                        y_initial,day_start_intervention,day_max);
                end
                I_intervention(j,:)=yi(:,2)';                
            end
            I_end_sample(:,k+1)=I_intervention(:,end);                       
            median_efficiency(i,(ni-1)*3+k)=median(I_intervention(:,end)./I_end_sample(:,1));
            opt_efficiency(i,(ni-1)*3+k)=I_intervention(ismember(idx,idx_opt),end)/I_end_sample(ismember(idx,idx_opt),1);            
           
            figure(f_I{ni});
            h=fill([t_all t_all(end:-1:1)]+delta_t,[min(I_intervention) max(I_intervention(:,end:-1:1))],...
                colors(k,:),'EdgeColor',[1 1 1]);
            set(h,'facealpha',.5)
            xlabel('Days since Feb 15, 2020');ylabel('Count');
            title(country_tags_unique{i});
        end
        scatter(48+delta_t,min(I_intervention(:,49))/2+max(I_intervention(:,49))/2,...
            80,colors(6,:),'p','filled','MarkerEdgeColor',[0 0 0]);
    end
end

% Visualize the intervention efficiency computed based on optimal parameter
% set (i.e. parameter set with optimal fitting to the data)
figure;
Intervention_names = {'All people wearing homemade masks',...
    'All people wearing surgical masks','All people wearing N95 respirators',...
    'Only infected individuals wearing homemade masks',...
    'Only infected individuals wearing surgical masks',...
    'Only infected individuals wearing N95 respirators'};
heatmap_cluster(opt_efficiency',Intervention_names,country_tags_unique,[0 1]);
title('Relative number of infected individuals on day 100 compared to no intervention');
colorbar;

%% Functions for simulation of different intervention strategies
function yi = SimIntervention(params_0,params_fc,y_initial,...
    day_start_intervention,day_max)
% Simulate disease dynamics under an intervention that changes parameters
% params_0 to params_0*params_fc
% t=0 to t=day_start_intervention: original parameters params_0
% t=day_start_intervention to t=day_max: altered parameters
% params_0*params_fc

% Simulate transmission dynamics before intervention
tstart=tic;
opts=odeset('Jacobian',@(t,y)modSIRD_Vpool_CBind_Jacobian(t,y,params_0),...
    'Events',@(t,y) Integration_Too_Slow(t,y,params_0,tstart),...
    'RelTol',1e-6,'AbsTol',1e-6);
[~,y_before]=ode23s(@(t,y)modSIRD_Vpool_CBind(t,y,params_0),...
    0:day_start_intervention,y_initial,opts);
y_transition=y_before(end,:); %System status on the day when intervention starts

% Simulate transmission dynamics after intervention
tstart=tic;
opts=odeset('Jacobian',@(t,y)modSIRD_Vpool_CBind_Jacobian(t,y,params_0.*params_fc),...
    'Events',@(t,y) Integration_Too_Slow(t,y,params_0,tstart),...
    'RelTol',1e-6,'AbsTol',1e-6);
[~,y_after]=ode23s(@(t,y)modSIRD_Vpool_CBind(t,y,params_0.*params_fc),...
    day_start_intervention:day_max,y_transition,opts);
yi=[y_before;y_after(2:end,:)];
end

function yo = SimDiseaseDynamics(params_0,y_initial,day_max)
% Simulate disease dynamics without intervention under parameters params_0
% Time range: t=0 to t=day_max

% Simulate transmission dynamics before intervention
tstart=tic;
opts=odeset('Jacobian',@(t,y)modSIRD_Vpool_CBind_Jacobian(t,y,params_0),...
    'Events',@(t,y) Integration_Too_Slow(t,y,params_0,tstart),...
    'RelTol',1e-6,'AbsTol',1e-6);
[~,yo]=ode23s(@(t,y)modSIRD_Vpool_CBind(t,y,params_0),...
    0:day_max,y_initial,opts);
end

%% Functions for the SIRD models
function dydt=modSIRD_Vpool_CBind(t,y,params)
%Reduced SIRD model with virus pool and cooperative infection
c=params(1);d=params(2);k=params(3);a=params(4);n=params(5);K=params(6);  %parameters
S=y(1);I=y(2);V=y(4); %variables
Hill_term=(V^n)/(V^n+K);
dydt=[-k*S*Hill_term;k*S*Hill_term-a*I;a*I;c*I-d*V];
end

function J=modSIRD_Vpool_CBind_Jacobian(t,y,params)
c=params(1);d=params(2);k=params(3);a=params(4);n=params(5);K=params(6);  %parameters
S=y(1);V=y(4); %variables
Hill_term=(V^n)/(V^n+K);
dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2))*n*V^(n-1);
if n==1
    dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2));
end
J=[-k*Hill_term 0 0 -k*S*dHdV;k*Hill_term -a 0 k*S*dHdV;0 a 0 0;0 c 0 -d];
end

function dydt=oriSIRD(t,y,params)
a=params(1);b=params(2);
% Reduced SIRD model
s=y(1);i=y(2); %variables
dydt=[-b*s*i;b*s*i-a*i;a*i];
end

function J=oriSIRD_Jacobian(t,y,params)
% Jacobian for the reduced SIRD model
a=params(1);b=params(2);
s=y(1);i=y(2); %variables
J=[-b*i -b*s 0;b*i b*s-a 0;0 a 0];
end

function dydt=modSIRD_Vpool(t,y,params)
% Reduced SIRD model with virus pool
c=params(1);d=params(2);k=params(3);a=params(4);  %parameters
S=y(1);I=y(2);V=y(4); %variables
dydt=[-k*S*V;k*S*V-a*I;a*I;c*I-d*V];
end

function J=modSIRD_Vpool_Jacobian(t,y,params)
% Jacobian matrix for reduced SIRD model with virus pool
c=params(1);d=params(2);k=params(3);a=params(4);  %parameters
S=y(1);V=y(4); %variables
J=[-k*V 0 0 -k*S;k*V -a 0 k*S;0 a 0 0;0 c 0 -d];
end

function [time_left,isterminal,direction] = Integration_Too_Slow(t,y,params,tstart)
% Event that happens if the integration has taken too much time. If this
% happens (integration takes longer than 1s), let the solver end prematurely
tc=toc(tstart);
time_left = max(1-tc,0);
isterminal=1;
direction=0;
if tc>=1
    t
    y
    params
    pause();
end
end
