%% Load results for state model fitting
daily_raw = readtable('daily.csv');
population_raw = readtable('population.csv');
states = intersect(unique(daily_raw.state),population_raw.Abbreviation);
params_opt_states = zeros(51,11);
fobj_opt_states = zeros(51,1);
label_good_fit = zeros(51,1);
params_samples_states = cell(51,1);
for i_state_sum=1:length(states)
    load(strcat('Fitting_results\',states{i_state_sum},'.mat'));
    day0 = find(total_cases_states(i_state_sum,:)~=0, 1);
    n_data_points = sum(~isnan(total_cases_states(i_state_sum,day0:end))) +...
        sum(~isnan(removed_states(i_state_sum,day0:end)));
    params_opt_states(i_state_sum,:) = x_opt';
    fobj_opt_states(i_state_sum) = f_opt/n_data_points;
    params_samples_states{i_state_sum} = x_kept;
    
    if size(x_kept,1)>1
        label_good_fit(i_state_sum)=1;
    end
end

%% Estimate death rate for each state
% Death rate is defined as the fraction of deaths in all removed cases
% Since most cases were severe on earlier dates, we assume that estimations
% based on later time points better reflect the actual values of death
% rates. Hence the last three dates with deaths and recoveries were used
% For states without data available for death rate estimation, average
% value of all other states was used
dr_dates = deceased_states./removed_states;
approx_death_rate_states = zeros(51,1);
for i=1:51
    if sum(isnan(dr_dates(i,:)))<size(dr_dates,2)
        x = dr_dates(i,:);
        x(isnan(x)) = [];
        approx_death_rate_states(i) = mean(x(end-2:end));
    end
end
approx_death_rate_states(approx_death_rate_states==0) = ...
    mean(approx_death_rate_states(approx_death_rate_states~=0));
clear dr_dates


%% Simulate the state models under three different schemes
n_days_sim = 501;
TimePoints = 1:10:501;
blues = brewermap(2,'Blues');

% Infection_Ratios: ratios of infections to population in states
% Death_Ratios: ratios of deaths to population in states
% Original: current intervention and reopen scenarios
% NoReopen: current intervention strength but no reopen
% NoIntervention: complete absence of intervention
Infection_Ratios_Original = -ones(51,length(TimePoints));
Death_Ratios_Original = -ones(51,length(TimePoints));
Infection_Ratios_NoReopen = -ones(51,length(TimePoints));
Infection_Ratios_NoIntervention = -ones(51,length(TimePoints));
Death_Ratios_NoReopen = -ones(51,length(TimePoints));
Death_Ratios_NoIntervention = -ones(51,length(TimePoints));

Infections_N_Original_All = zeros(500,n_days_sim,51);
%Numbers of infections simulated for all states, each with 500 parameter
%sets
Deaths_N_Original_All = zeros(500,n_days_sim,51);
%Numbers of deaths simulated for all states, each with 500 parameter sets
Infections_N_NoReopen_All = zeros(500,n_days_sim,51);
Deaths_N_NoReopen_All = zeros(500,n_days_sim,51);
Infections_N_NoIntervention_All = zeros(500,n_days_sim,51);
Deaths_N_NoIntervention_All = zeros(500,n_days_sim,51);

for i_state = 1:51
    [~,idx_population] = ismember(states{i_state},population_raw.Abbreviation);
    day0 = find(total_cases_states(i_state,:)~=0, 1 );
    state_names(i_state)
    if label_good_fit(i_state) == 1
        params_original_samp = params_samples_states{i_state};
        params_original_samp = [params_original_samp(randperm(size(params_original_samp,1),499),:);...
            params_opt_states(i_state,:)];
        params_noreopen_samp = params_original_samp;
        params_noreopen_samp(:,10) = 4;
        params_nointervention_samp = params_original_samp;
        params_nointervention_samp(:,9) = 4;
        List_params_samp_scenarios = {params_original_samp,params_noreopen_samp,params_nointervention_samp};
        %Simulate the original model (i.e. with the current scenario of
        %non-pharmaceutical intervention
        %Add data for dates before identification of the 1st case in that
        %state
        for i_scenario = 1:3
            params_samp = List_params_samp_scenarios{i_scenario};
            I_time_course = zeros(500,n_days_sim);
            N_time_course = zeros(500,n_days_sim);
            for i = 1:500
                time_course = [zeros(day0-1,4);simSIRDModel(@cooperative_SIR_DI,@cooperative_SIR_DI_Jacobian,...
                    n_days_sim+1-day0,params_samp(i,:),total_cases_states(i_state,:),...
                    removed_states(i_state,:),population_states(idx_population))];
                I_time_course(i,:) = time_course(:,2)';
                N_time_course(i,:) = time_course(:,3)';
            end
            switch i_scenario
                case 1
                    Infections_N_Original_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_Original_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_Original(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_Original(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
                    
                    % Save the simulation results to figure
                    figure;
                    subplot(1,2,1);hold on;
                    fill([1:n_days_sim n_days_sim:-1:1],...
                        [min(I_time_course+N_time_course) max(I_time_course(:,end:-1:1)...
                        +N_time_course(:,end:-1:1))],blues(1,:),'EdgeColor',[1 1 1]);
                    plot(1:n_days_sim,I_time_course(end,:)+N_time_course(end,:),'Color',blues(2,:));
                    scatter(1:size(total_cases_states,2),total_cases_states(i_state,:),100,blues(2,:),'.');
                    title(sprintf('%s %s','Infections in',state_names{i_state}));
                    xlim([1 501]);xticks(1:50:501);
                    xticklabels(date2str(1:50:501));xtickangle(45);
                    ylabel('Number of infections');
                    
                    subplot(1,2,2);hold on;
                    fill([1:n_days_sim n_days_sim:-1:1],[min(N_time_course) max(N_time_course(:,end:-1:1))],...
                        blues(1,:),'EdgeColor',[1 1 1]);
                    plot(1:n_days_sim,N_time_course(end,:),'Color',blues(2,:));
                    scatter(1:size(total_cases_states,2),removed_states(i_state,:),100,blues(2,:),'.');
                    title(sprintf('%s %s','Removed cases in',state_names{i_state}));
                    xlim([1 501]);xticks(1:50:501);
                    xticklabels(date2str(1:50:501));xtickangle(45);
                    ylabel('Number of recovered/deceased people');
                    
                    saveas(gcf,strcat(state_names{i_state},'.fig'));
                case 2
                    Infections_N_NoReopen_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_NoReopen_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_NoReopen(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_NoReopen(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
                case 3
                    Infections_N_NoIntervention_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_NoIntervention_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_NoIntervention(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_NoIntervention(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
            end
        end
    else
        params_original_samp = params_opt_states(i_state,:);
        params_noreopen_samp = params_original_samp;
        params_noreopen_samp(10) = 4;
        params_nointervention_samp = params_original_samp;
        params_nointervention_samp(9) = 4;
        List_params_samp_scenarios = {params_original_samp,params_noreopen_samp,params_nointervention_samp};
        
        for i_scenario = 1:3
            params_samp = List_params_samp_scenarios{i_scenario};
            time_course = [zeros(day0-1,4);simSIRDModel(@cooperative_SIR_DI,@cooperative_SIR_DI_Jacobian,...
                n_days_sim+1-day0,params_samp,total_cases_states(i_state,:),...
                removed_states(i_state,:),population_states(idx_population))];
            I_time_course = repmat(time_course(:,2)',500,1);
            N_time_course = repmat(time_course(:,3)',500,1);
            switch i_scenario
                case 1
                    Infections_N_Original_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_Original_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_Original(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_Original(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
                    
                    % Save the simulation results to figure
                    figure;
                    subplot(1,2,1);hold on;
                    fill([1:n_days_sim n_days_sim:-1:1],...
                        [min(I_time_course+N_time_course) max(I_time_course(:,end:-1:1)...
                        +N_time_course(:,end:-1:1))],blues(1,:),'EdgeColor',[1 1 1]);
                    plot(1:n_days_sim,I_time_course(end,:)+N_time_course(end,:),'Color',blues(2,:));
                    scatter(1:size(total_cases_states,2),total_cases_states(i_state,:),100,blues(2,:),'.');
                    title(sprintf('%s %s','Infections in',state_names{i_state}));
                    xlim([1 501]);xticks(1:50:501);
                    xticklabels(date2str(1:50:501));xtickangle(45);
                    ylabel('Number of infections');
                    
                    subplot(1,2,2);hold on;
                    fill([1:n_days_sim n_days_sim:-1:1],[min(N_time_course) max(N_time_course(:,end:-1:1))],...
                        blues(1,:),'EdgeColor',[1 1 1]);
                    plot(1:n_days_sim,N_time_course(end,:),'Color',blues(2,:));
                    scatter(1:size(total_cases_states,2),removed_states(i_state,:),100,blues(2,:),'.');
                    title(sprintf('%s %s','Removed cases in',state_names{i_state}));
                    xlim([1 501]);xticks(1:50:501);
                    xticklabels(date2str(1:50:501));xtickangle(45);
                    ylabel('Number of recovered/deceased people');       
                    saveas(gcf,strcat(state_names{i_state},'.fig'));
                case 2
                    Infections_N_NoReopen_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_NoReopen_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_NoReopen(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_NoReopen(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
                case 3
                    Infections_N_NoIntervention_All(:,:,i_state) = ...
                        I_time_course + N_time_course;
                    Deaths_N_NoIntervention_All(:,:,i_state) = ...
                        N_time_course*approx_death_rate_states(i_state);
                    for k=1:length(TimePoints)
                        Infection_Ratios_NoIntervention(i_state,k) = (I_time_course(end,TimePoints(k))+...
                            N_time_course(end,TimePoints(k)))/population_states(idx_population);
                        Death_Ratios_NoIntervention(i_state,k) = N_time_course(end,TimePoints(k))...
                            *approx_death_rate_states(i_state)/population_states(idx_population);
                    end
            end
        end
    end
    
    %{
    
    %}
end

Infections_reduced_intervention = Infection_Ratios_NoIntervention - Infection_Ratios_Original;
Infections_increased_reopen = Infection_Ratios_Original - Infection_Ratios_NoReopen;
Lives_saved_intervention = Death_Ratios_NoIntervention - Death_Ratios_Original;
Lives_cost_reopen = Death_Ratios_Original - Death_Ratios_NoReopen;

Infections_reduced_intervention(Infections_reduced_intervention<0) = 0;
Infections_increased_reopen(Infections_increased_reopen<0) = 0;
Lives_saved_intervention(Lives_saved_intervention<0) = 0;
Lives_cost_reopen(Lives_cost_reopen<0) = 0;

clear i_scenario i_state i_state_sum params_samp I_time_course N_time_course k

%% Visualize results on the USA map
titles = date2str(TimePoints);

Data = {Infection_Ratios_Original,Infection_Ratios_NoReopen,Infection_Ratios_NoIntervention,...
    Death_Ratios_Original,Death_Ratios_NoReopen,Death_Ratios_NoIntervention};
RunNames = {'Infections_Current','Infections_NoReopen','Infections_NoIntervention',...
    'Deaths_Current','Deaths_NoReopen','Deaths_NoIntervention'};
var_limits = {[0 1],[0 1],[0 1],[0 0.1],[0 0.1],[0 0.1]};

for i=1:6
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf, 'color', 'white');
    Draw_on_map(state_names,Data{i},RunNames{i},titles,var_limits{i});
end

%% Plot intervention strength values on US map
figure;
set(gcf, 'color', 'white');
Draw_on_map_static(state_names,min_alpha_states,[0 1]);
figure;
set(gcf, 'color', 'white');
Draw_on_map_static(state_names,mean_alpha_states,[0 1]);

%% Calculate absolute numbers of infections/deaths reduced by the
% intervention or increased by reopen
Total_Infections_Current_Sample = sum(Infections_N_Original_All,3);
Total_Infections_NoReopen_Sample = sum(Infections_N_NoReopen_All,3);
Total_Infections_NoIntervention_Sample = sum(Infections_N_NoIntervention_All,3);
Total_Deaths_Current_Sample = sum(Deaths_N_Original_All,3);
Total_Deaths_NoReopen_Sample = sum(Deaths_N_NoReopen_All,3);
Total_Deaths_NoIntervention_Sample = sum(Deaths_N_NoIntervention_All,3);
% Create lists for iteration
List_Infections = {Total_Infections_Current_Sample,Total_Infections_NoIntervention_Sample,...
    Total_Infections_NoReopen_Sample};
List_Deaths = {Total_Deaths_Current_Sample,Total_Deaths_NoIntervention_Sample,...
    Total_Deaths_NoReopen_Sample};

% Define colors
colors = brewermap(8,'Set2');

figure;
subplot(1,2,1);
hold on;
for i_list = 1:3
    data = List_Infections{i_list};
    h=fill([0:500 500:-1:0],[min(data) max(data(:,end:-1:1))],colors(i_list,:),'EdgeColor',[1 1 1]);
    set(h,'facealpha',.5)
    plot(0:500,data(end,:),'Color',colors(i_list,:));
end
xticks(0:100:500);xticklabels(date2str(1:100:501));
xtickangle(45);xlim([0 500]);
ylabel('Number of people');title('Infections');
legend('Current intervention scenarios (500 samples)','Current intervention scenarios (optimal fitting)',...
    'No intervention (500 samples)','No intervention (optimal fitting)',...
    'No reopen (500 samples)','No reopen (optimal fitting)');
subplot(1,2,2);
hold on;
for i_list = 1:3
    data = List_Deaths{i_list};
    h=fill([0:500 500:-1:0],[min(data) max(data(:,end:-1:1))],colors(i_list,:),'EdgeColor',[1 1 1]);
    set(h,'facealpha',.5)
    plot(0:500,data(end,:),'Color',colors(i_list,:));
end
xticks(0:100:500);xticklabels(date2str(1:100:501));
xtickangle(45);xlim([0 500]);
ylabel('Number of people');title('Deaths');
legend('Current intervention scenarios','Current intervention scenarios',...
    'No intervention','No intervention',...
    'No reopen','No reopen');
clear i_list data

%% Compare intervention strength across states
t = 1:size(total_cases_states,2); % Examine the time window from Jan 22 to July 15, 2020
alpha = zeros(51,size(total_cases_states,2));
for i=1:51
    day0 = find(total_cases_states(i,:)>0, 1);
    alpha(i,:) = arrayfun(@(x)intervention_strength(x-day0,10.^params_opt_states(i,:)),t);
end
mean_alpha_states = mean(alpha,2);
min_alpha_states = min(alpha,[],2);
x=[mean_alpha_states min_alpha_states];
figure;violinplot(x,{'Average','Minimal value'});
title('Average and minimal \alpha values in U.S. states')
ylabel('\alpha');box on;

%% Functions for model definition and result visualization
function Draw_on_map_static(states,values,val_limits)
%Visualize data in values (corresponding to states in the input variable
%'states') on US map
ax = usamap({'WA','FL','ME','AK'});
set(ax, 'Visible', 'off')
latlim = getm(ax, 'MapLatLimit');
lonlim = getm(ax, 'MapLonLimit');

struct_states_map = shaperead('usastatehi','UseGeoCoords', true);
state_abbreviations = {'AL','AK','AZ','AR','CA','CO','CT','DE','FL','GA',...
    'HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS',...
    'MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA',...
    'RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY','DC'};
state_names_map = {struct_states_map.Name};
[~,idx] = ismember(state_names_map,states);
for i=1:51
    struct_states_map(i).ShowValue = 0;
    struct_states_map(i).Abbreviation = state_abbreviations{i};
end
for i=1:length(idx)
    struct_states_map(i).ShowValue = values(idx(i));
end
cmap = brewermap(64,'RdYlBu');
cmap = cmap(end:-1:1,:);
%val_limits = [min(values) max(values)];
%val_limits = [0 1];

colormap(cmap);
caxis(val_limits);
state_colors = interp1(linspace(val_limits(1),val_limits(2),size(cmap,1)),...
    cmap,values(idx));
surfaceColors = makesymbolspec('Polygon',{'INDEX',[1 51],...
    'FaceColor',state_colors});
geoshow(struct_states_map,'DisplayType','polygon','SymbolSpec',surfaceColors);
colorbar;
%{
% Code for displaying state abbreviations
lat = [struct_states_map.LabelLat];
lon = [struct_states_map.LabelLon];
tf = ingeoquad(lat, lon, latlim, lonlim);
textm(lat(tf), lon(tf), {struct_states_map(tf).Abbreviation}, ...
    'HorizontalAlignment', 'center')
%}
end

function s = date2str(date_idx)
% Generate textual identifiers for dates
% dates start from Jan 22, 2020
day_start = datetime(2020,1,22);
day_now = dateshift(day_start,'end','day',date_idx-2);
s = string(datestr(day_now));
end

function Draw_on_map(states,values,draw_name,draw_titles,val_limits)
%Visualize data in values (corresponding to states in the input variable
%'states') on US map
ax = usamap({'WA','FL','ME','AK'});
set(ax, 'Visible', 'off')
latlim = getm(ax, 'MapLatLimit');
lonlim = getm(ax, 'MapLonLimit');

struct_states_map = shaperead('usastatehi','UseGeoCoords', true);
state_abbreviations = {'AL','AK','AZ','AR','CA','CO','CT','DE','FL','GA',...
    'HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS',...
    'MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA',...
    'RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY','DC'};
state_names_map = {struct_states_map.Name};
[~,idx] = ismember(state_names_map,states);
for i=1:51
    struct_states_map(i).ShowValue = 0;
    struct_states_map(i).Abbreviation = state_abbreviations{i};
end
for i=1:length(idx)
    struct_states_map(i).ShowValue = values(idx(i));
end
cmap = brewermap(64,'RdYlBu');
cmap = cmap(end:-1:1,:);
%val_limits = [min(values) max(values)];
%val_limits = [0 1];

if nargin == 4
    val_limits = [0 max(values(:))];
end
colormap(cmap);
caxis(val_limits);

gif_name = strcat(draw_name,'.gif');
for j=1:size(values,2)
    state_colors = interp1(linspace(val_limits(1),val_limits(2),size(cmap,1)),...
        cmap,values(idx,j));
    state_colors(values(idx,j)==-1,:)=repmat([1 1 1],sum(values(:,j)==-1),1);
    
    surfaceColors = makesymbolspec('Polygon',{'INDEX',[1 51],...
        'FaceColor',state_colors});
    
    geoshow(struct_states_map,'DisplayType','polygon','SymbolSpec',surfaceColors);
    
    colorbar;
    %{
    lat = [struct_states_map.LabelLat];
    lon = [struct_states_map.LabelLon];
    tf = ingeoquad(lat, lon, latlim, lonlim);
    textm(lat(tf), lon(tf), {struct_states_map(tf).Abbreviation}, ...
        'HorizontalAlignment', 'center')
    %}
    title(draw_titles{j});
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if j == 1
        imwrite(imind,cm,gif_name,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,gif_name,'gif','WriteMode','append');
    end
end
end

function alpha = intervention_strength(t,params)
%Compute time-dependent intervention strength alpha
A=params(7);B=params(8);T0=params(9);Delta=params(10);N=params(11);
if t>T0
    t_rel = t-T0;
    t_rel_max = N/B;
    if t_rel<t_rel_max
        alpha = 1/(1+A*(t_rel^N)*exp(-B*t_rel));
    elseif t_rel<t_rel_max+Delta
        alpha = 1/(1+A*((N/B)^N)*exp(-N));
    else
        alpha = 1/(1+A*((t_rel-Delta)^N)*exp(-B*(t_rel-Delta)));
    end
else
    alpha = 1;
end
end

function dydt=cooperative_SIR_DI(t,y,params)
%SIRD model with virus pool and cooperativity in infection by virus,
%considering dynamic intervention (i.e. intervention strength alpha changes
%with time t by the function alpha = 1/(1+A*t*exp(-B*t))
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);K=params(6);  %parameters
A=params(7);B=params(8);T0=params(9);Delta=params(10);N=params(11);

if t>T0
    t_rel = t-T0;
    t_rel_max = N/B;
    if t_rel<t_rel_max
        alpha = 1/(1+A*(t_rel^N)*exp(-B*t_rel));
    elseif t_rel<t_rel_max+Delta
        alpha = 1/(1+A*((N/B)^N)*exp(-N));
    else
        alpha = 1/(1+A*((t_rel-Delta)^N)*exp(-B*(t_rel-Delta)));
    end
    c=c*alpha;
    K=K/(alpha^n);
end

S=y(1);I=y(2);V=y(4); %variables
Hill_term=(V^n)/(V^n+K);
dydt=[-k*S*Hill_term;k*S*Hill_term-lambda*I;lambda*I;c*I-d*V];
end

function J=cooperative_SIR_DI_Jacobian(t,y,params)
% Jacobian for cooperative SIR model, considering dynamic intervention
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);K=params(6);  %parameters
S=y(1);V=y(4); %variables

A=params(7);B=params(8);T0=params(9);Delta=params(10);N=params(11);
if t>T0
    t_rel = t-T0;
    t_rel_max = N/B;
    if t_rel<t_rel_max
        alpha = 1/(1+A*(t_rel^N)*exp(-B*t_rel));
    elseif t_rel<t_rel_max+Delta
        alpha = 1/(1+A*((N/B)^N)*exp(-N));
    else
        alpha = 1/(1+A*((t_rel-Delta)^N)*exp(-B*(t_rel-Delta)));
    end
    c=c*alpha;
    K=K/(alpha^n);
end

Hill_term=(V^n)/(V^n+K);
dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2))*n*V^(n-1);
if n==1
    dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2));
end
J=[-k*Hill_term 0 0 -k*S*dHdV;k*Hill_term -lambda 0 k*S*dHdV;0 lambda 0 0;0 c 0 -d];
end

function [time_left,isterminal,direction] = Integration_Too_Slow(t,y,params,tstart)
% Event that happens if the integration has taken too much time. If this
% happens (integration takes longer than 1s), let the solver end prematurely
tc=toc(tstart);
time_left = max(1-tc,0);
isterminal=1;
direction=0;
end

%% Function for model simulation
function time_course=simSIRDModel(Func_ODE,Func_Jacobian,...
    n_days_simulation,log_params,data_cases_total,data_cases_removed,...
    total_population)
%Function for fitting a specific pandemic model (dydtModelFunc,
%JacobianModelFunc) to population data
% Input variables:
% - 1. n_days_simulation: number of days to be simulated
% - 2. Func_ODE: function handle for the ODEs model
% - 3. Func_Jacobian: function handle for the Jacobian of the ODE model
% - 4. params: parameters used for simulation
% - 5. data_cases_total: known numbers of total cases
% - 6. data_cases_removed: known numbers of removed cases
% - 7. total_population: total population in the state

day0 = find(data_cases_total~=0, 1 ); % Find the first day with reported case
days_prediction = day0:day0+n_days_simulation-1;
data_cases_active = data_cases_total - data_cases_removed;
if isnan(data_cases_active(day0))
    y_initial = [total_population - data_cases_total(day0) data_cases_total(day0) 0 0];
else
    y_initial=[total_population-data_cases_total(day0) data_cases_active(day0) data_cases_removed(day0) 0];
end
params = 10.^log_params;

tstart=tic;
opts=odeset('Jacobian',@(t,y)Func_Jacobian(t,y,params),...
    'Events',@(t,y) Integration_Too_Slow(t,y,params,tstart),'RelTol',1e-6,'AbsTol',1e-6);
[~,time_course]=ode23s(@(t,y)Func_ODE(t,y,params),days_prediction - day0,y_initial,opts);
time_course = real(time_course);
end
