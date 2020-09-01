% Read daily COVID-19 case data in US states
daily_raw = readtable('daily.csv');
population_raw = readtable('population.csv');
states = intersect(unique(daily_raw.state),population_raw.Abbreviation);
dates = unique(daily_raw.date);
[~,idx_states]=ismember(states,population_raw.Abbreviation);
state_names = population_raw.State(idx_states);

[~,idx_states]=ismember(daily_raw.state,states);
[~,idx_dates]=ismember(daily_raw.date,dates);
idx_keep = find(idx_states>0);

% Generate state-by-state data matrice for total, deceased and recovered
% SARS-nCov-2 case counts
total_cases_states = full(sparse(idx_states(idx_keep),idx_dates(idx_keep),...
    daily_raw.positive(idx_keep)));
deceased_states = full(sparse(idx_states(idx_keep),idx_dates(idx_keep),...
    daily_raw.death(idx_keep)));
recovered_states = full(sparse(idx_states(idx_keep),idx_dates(idx_keep),...
    daily_raw.recovered(idx_keep)));

% Preprocess recovered_states to correct for bias caused by weekly reports
% (i.e. unchanged numbers of recovered cases during the time window of one
% week)
for i=1:length(states)
    j_start = find(recovered_states(i,:)>0, 1);
    while j_start < length(dates)
        j_end = find(recovered_states(i,:) == recovered_states(i,j_start),1,'last');
        if j_end > j_start
            recovered_states(i,j_start+1:j_end) = NaN;
            j_start = j_end + 1;
        else
            j_start = j_start + 1;
        end
    end
end

removed_states = deceased_states + recovered_states;

clear idx_keep idx_states idx_dates

%% Select states with most serious covid-19 spread recently
% Define state-by-state metrics for covid-19 spread
%{
latest_daily_increase = total_cases_states(:,end)-total_cases_states(:,end-1);
[~,rank_states_by_increase] = sort(latest_daily_increase,'descend');
n_keep = 10;
figure;
subplot(1,2,1);
heatmap_cluster_1d(total_cases_states(rank_states_by_increase(1:n_keep),:)./...
    total_cases_states(rank_states_by_increase(1:n_keep),end),...
    state_names(rank_states_by_increase(1:n_keep)),[],[0 1],1);
subplot(1,2,2);
plot(total_cases_states(rank_states_by_increase(1:n_keep),:)');

figure;
Draw_on_map(state_names,latest_daily_increase);
%}
%% Fit models to state data
%Use the model that considers under-reporting
params_lb = [1e-8 1 1e-8 1e-6 1.1 1e-3 1e-4 1e-4 1 1 0.1];
params_ub = [1e-2 100 1e-2 1e2 5 1e2 100 100 100 1000 2];

%i_state = rank_states_by_increase(2);
[~,idx_population] = ismember(states{i_state},population_raw.Abbreviation);
population_states = population_raw.Population;
states{i_state}
[x_opt,f_opt,x_kept,f_kept]=fitSIRDModel(@cooperative_SIR_DI,@cooperative_SIR_DI_Jacobian,...
    total_cases_states(i_state,:),removed_states(i_state,:),...
    population_states(idx_population),params_lb,params_ub,states{i_state});
save(strcat(states{i_state},'.mat'));
exit;

%% Functions for model definition and result visualization
function Draw_on_map(states,values)
%Visualize data in values (corresponding to states in the input variable
%'states') on US map
ax = usamap({'WA','FL','AK','ME'});
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
%{
val_limits = [min(values) max(values)];
surfaceColors = makesymbolspec('Polygon',{'INDEX',[1 51],...
    'FaceColor',interp1(linspace(val_limits(1),val_limits(2),size(cmap,1)),...
    cmap,values(idx))});
%}
[~,value_orders] = sort(values);
[~,value_ranks] = sort(value_orders);
surfaceColors = makesymbolspec('Polygon',{'INDEX',[1 51],...
    'FaceColor',interp1(linspace(1,51,size(cmap,1)),...
    cmap,value_ranks(idx))});

geoshow(struct_states_map,'DisplayType','polygon','SymbolSpec',surfaceColors);

lat = [struct_states_map.LabelLat];
lon = [struct_states_map.LabelLon];
tf = ingeoquad(lat, lon, latlim, lonlim);
textm(lat(tf), lon(tf), {struct_states_map(tf).Abbreviation}, ...
    'HorizontalAlignment', 'center')

end

function dydt=cooperative_SIR(t,y,params)
%SIRD model with virus pool and cooperativity in infection by virus
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);K=params(6);  %parameters
S=y(1);I=y(2);V=y(4); %variables
Hill_term=(V^n)/(V^n+K);
dydt=[-k*S*Hill_term;k*S*Hill_term-lambda*I;lambda*I;c*I-d*V];
end

function J=cooperative_SIR_Jacobian(t,y,params)
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);K=params(6);  %parameters
S=y(1);V=y(4); %variables
Hill_term=(V^n)/(V^n+K);
dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2))*n*V^(n-1);
if n==1
    dHdV=(1/(V^n+K)-(V^n)/((V^n+K)^2));
end
J=[-k*Hill_term 0 0 -k*S*dHdV;k*Hill_term -lambda 0 k*S*dHdV;0 lambda 0 0;0 c 0 -d];
end

%%
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
%{
if tc>=1
    t
    y
    params
    pause();
end
%}
end

%% Function for model fitting
function [x_opt,f_opt,x_kept,f_kept]=fitSIRDModel(Func_ODE,Func_Jacobian,...
    data_cases_total,data_cases_removed,...
    total_population,params_lb,params_ub,model_name)
%Function for fitting a specific pandemic model (dydtModelFunc,
%JacobianModelFunc) to population data
% Input variables:
% - 1. data_cases_total: time course data of total infections
% - 2. data_cases_removed: time course data for removed (recovered or
% deceased) cases
% - 3. total_population: total population in the state
% - 4. params_lb: lower bounds for parameters in the model
% - 5. params_ub: upper bounds for parameters in the model

day0 = find(data_cases_total~=0, 1 ); % Find the first day with reported case
day_max = length(data_cases_total);%170; % Last date used for model fitting
days_data_available = day0:day_max;
days_data_validation = day_max+1:length(data_cases_total);
days_prediction = day0:length(data_cases_total);
data_cases_active = data_cases_total - data_cases_removed;

if isnan(data_cases_active(day0))
    y_initial = [total_population - data_cases_total(day0) data_cases_total(day0) 0 0];
else
    y_initial=[total_population-data_cases_total(day0) data_cases_active(day0) data_cases_removed(day0) 0];
end

%Retrieve data points for fitting and store them in a 2xn matrix in which
%the first row stores total case numbers and second row stores numbers of
%removed cases
data_fitting = [data_cases_total(days_data_available);data_cases_removed(days_data_available)];
data_validation = [data_cases_total(days_data_validation);data_cases_removed(days_data_validation)];
n_data_points = sum(~isnan(data_fitting(:))); % Number of data points available for fitting

% Fit the cooperative infection model to data
DSA_options.start_beta=1e-6;
DSA_options.end_beta=1e8;
DSA_options.expected_fobj=n_data_points*0.00001;%0.00005;
DSA_options.visual_output=0;

% Options used for debugging in which running time will be greatly reduced
% by decreasing Markov chain length and increasing cooling rate
%{
DSA_options.cooling_rate = 1.1;
DSA_options.mc_length = 100;
%}

[x_opt,f_opt,x_kept,f_kept] = DSA(@calResidues,log10(params_lb(:)),...
    log10(params_ub(:)),DSA_options);

I_samp=zeros(size(x_kept,1),length(days_prediction)); %Total infections (I+N in the model)
N_samp=zeros(size(x_kept,1),length(days_prediction)); %Removed cases (N in the model)
for i=1:size(x_kept,1)
    params=10.^x_kept(i,:);
    opts=odeset('Jacobian',@(t,y)Func_Jacobian(t,y,params));
    [t,y]=ode23s(@(t,y)Func_ODE(t,y,params),days_prediction-day0,y_initial,opts);
%{    
    if length(params_lb)==100 %Model with consideration of underreporting
        a = params(7);
        b = params(8);
        c = params(9);
        IN_actual = [0;y(:,2)+y(:,3)];
        N_actual = [0;y(:,3)];
        p_report = 1./(a+b*exp(-c*(days_prediction-day0)));
        dIN_report = (IN_actual(2:end)-IN_actual(1:end-1)).*p_report(:);
        dN_report = (N_actual(2:end)-N_actual(1:end-1)).*p_report(:);
        IN_report = cumsum(dIN_report);
        N_report = cumsum(dN_report);
        I_samp(i,:) = IN_report';
        N_samp(i,:) = N_report';
    else
    %}
        I_samp(i,:)=y(:,2)'+y(:,3)';
        N_samp(i,:)=y(:,3)';
    %end
end
I_lb=min(I_samp);I_ub=max(I_samp);
N_lb=min(N_samp);N_ub=max(N_samp);

params_opt=10.^x_opt;
opts=odeset('Jacobian',@(t,y)Func_Jacobian(t,y,params_opt));
[t,y]=ode23s(@(t,y)Func_ODE(t,y,params_opt),days_prediction-day0,y_initial,opts);
%{
if length(params_lb)==100 %Model with consideration of underreporting
    a = params_opt(7);
    b = params_opt(8);
    c = params_opt(9);
    IN_actual = [0;y(:,2)+y(:,3)];
    N_actual = [0;y(:,3)];
    p_report = 1./(a+b*exp(-c*(days_prediction-day0)));
    dIN_report = (IN_actual(2:end)-IN_actual(1:end-1)).*p_report(:);
    dN_report = (N_actual(2:end)-N_actual(1:end-1)).*p_report(:);
    IN_report = cumsum(dIN_report);
    N_report = cumsum(dN_report);
    y(:,2) = IN_report - N_report;
    y(:,3) = N_report;
end
%}
if size(x_kept,1)==0
    I_lb=y(:,2)'+y(:,3)';
    I_ub=y(:,2)'+y(:,3)';
    N_lb=y(:,3)';
    N_ub=y(:,3)';
end

figure;
subplot(1,2,1);
fill([t' t(end:-1:1)'],[I_lb I_ub(end:-1:1)],[0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on;plot(t,y(:,2)+y(:,3));scatter(days_data_available-day0,data_fitting(1,:),'o');
scatter(days_data_validation-day0,data_validation(1,:),'o');
title(strcat(model_name,':Total infected'));
xlabel('Days');ylabel('Count');
subplot(1,2,2);
fill([t' t(end:-1:1)'],[N_lb N_ub(end:-1:1)],[0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on;plot(t,y(:,3));scatter(days_data_available-day0,data_fitting(2,:),'o');
scatter(days_data_validation-day0,data_validation(2,:),'o');
title(strcat(model_name,':Recovered/deceased'));
xlabel('Days');ylabel('Count');
saveas(gcf,strcat(model_name,'.fig'));

    function f=calResidues(log_params)
        params=10.^log_params;
        tstart=tic;
        %max(data_fitting')
        
        opts=odeset('Jacobian',@(t,y)Func_Jacobian(t,y,params),...
            'Events',@(t,y) Integration_Too_Slow(t,y,params,tstart),'RelTol',1e-6,'AbsTol',1e-6);
        [t,y]=ode23s(@(t,y)Func_ODE(t,y,params),days_data_available-day0,y_initial,opts);
        
        if size(y,1) ~= length(days_data_available)
            f=10000*ones(n_data_points,1);
        else
           % if 1%length(log_params)~=100 %Model without consideration of underreporting
                f = ([y(:,2)+y(:,3) y(:,3)] - data_fitting')./max(data_fitting');
                f = f(:);                             
                f(isnan(f)) = [];
                %{
            else %Model with consideration of underreporting
                a = params(7);
                b = params(8);
                c = params(9);
                IN_actual = [0;y(:,2)+y(:,3)];
                N_actual = [0;y(:,3)];
                p_report = 1./(a+b*exp(-c*(days_data_available-day0)));
                dIN_report = (IN_actual(2:end)-IN_actual(1:end-1)).*p_report(:);
                dN_report = (N_actual(2:end)-N_actual(1:end-1)).*p_report(:);
                IN_report = cumsum(dIN_report);
                N_report = cumsum(dN_report);
                f = ([IN_report N_report]-data_fitting')./max(data_fitting');
                f = f(:);
                f(isnan(f))=[];
            end
                %}
            f=real(f);
        end
        
        %For troubleshooting
        if norm(f)==0
            p_report
            [a b c]
            params
            pause();
        end
    end
end
