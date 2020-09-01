%--------------------------------------------------------------------------------------
% Script for fitting country-specific COVID models
%--------------------------------------------------------------------------------------
% Arguments needed:
%
% i_country - index of the country for which the parameter estimation is done
% Values of i_country:
% 1 - USA
% 2 - Italy
% 3 - Spain
% 4 - Germany
% 5 - France
% 6 - Iran
% 7 - UK
% 8 - Turkey
% 9 - Switzerland
%
% i_model - index of the model used in fitting the transmission data
% Values of i_model:
% 1 - reduced SIRD models
% 2 - reduced SIRD model with virus pool
% 3 - reduced SIRD model with virus pool and cooperative infection
%--------------------------------------------------------------------------------------
% Input files:
% TotalCases.csv - time-course data for total case numbers in the nine countries
% ActiveCases.csv - time-course data for active case numbers in the nine countries
%--------------------------------------------------------------------------------------

% Add the folder 'Functions' to the MATLAB search path
addpath('Functions');

% Define parameters and read transmission dynamics data
TotalPopulations = [329065000 60550000 46737000 83517000 65130000 ...
    82914000 67530000 83430000 8591000]; %Population (2019 UN data) in the nine countries
ActiveCases = readtable('ActiveCases.csv','ReadRowNames',true);
TotalCases = readtable('TotalCases.csv','ReadRowNames',true);
Countries = TotalCases.Properties.RowNames;
Dates = TotalCases.Properties.VariableNames;

% Identify the onset of the pandemic
if ActiveCases{i_country,1}>0
    start_idx = 1;
else
    start_idx = max(find(ActiveCases{i_country,:}==0))+1;
end
start_date = Dates{start_idx};

Models = {'SIRD original','SIRD virus pool','SIRD cooperative infection'}; % Names of models
% Lower and upper bounds for model parameters
Params_LBs = {[1e-6 1e-12],[1e-8 1 1e-8 1e-6],...
    [1e-8 1 1e-8 1e-6 1.1 1e-3]};
Params_UBs = {[1e2 1e-6],[1e-2 100 1e-2 1e2],...
    [1e-2 100 1e-2 1e2 5 1e2]};
ODE_Funcs = {@oriSIRD,@modSIRD_Vpool,@modSIRD_Vpool_CBind}; % Function handles for the models
Jacobian_Funcs = {@oriSIRD_Jacobian,@modSIRD_Vpool_Jacobian,@modSIRD_Vpool_CBind_Jacobian}; % Function handles for Jacobian of the models
Virus_Tags = [0 1 1]; % Binary variable indicating whether the model includes virus term
Days_Prediction = 0:80; % Time span for the simulation

run_name = strcat(Countries{i_country},'-',Models{i_model});
[x_opt,f_opt,x_kept,f_kept]=fitSIRDModel(ODE_Funcs{i_model},Jacobian_Funcs{i_model},...
run_name,Virus_Tags(i_model),ActiveCases{i_country,start_idx:end},...
TotalCases{i_country,start_idx:end},TotalPopulations(i_country),...
    Params_LBs{i_model},Params_UBs{i_model},Days_Prediction,start_date);
save(strcat(run_name,'.mat'));
exit;

%% Definition of functions
function dydt=modSIRD_Vpool_CBind(t,y,params)
%SIRD model with virus pool and cooperativity in infection by virus
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);b=params(6);  %parameters
S=y(1);I=y(2);V=y(4); %variables
Hill_term=(V^n)/(V^n+b);
dydt=[-k*S*Hill_term;k*S*Hill_term-lambda*I;lambda*I;c*I-d*V];
end

function J=modSIRD_Vpool_CBind_Jacobian(t,y,params)
c=params(1);d=params(2);k=params(3);lambda=params(4);n=params(5);b=params(6);  %parameters
S=y(1);V=y(4); %variables
Hill_term=(V^n)/(V^n+b);
dHdV=(1/(V^n+b)-(V^n)/((V^n+b)^2))*n*V^(n-1);
if n==1
    dHdV=(1/(V^n+b)-(V^n)/((V^n+b)^2));
end
J=[-k*Hill_term 0 0 -k*S*dHdV;k*Hill_term -lambda 0 k*S*dHdV;0 lambda 0 0;0 c 0 -d];
end

function dydt=oriSIRD(t,y,params)
a=params(1);b=params(2);
% ODEs for the original SIRD model
s=y(1);i=y(2);n=y(3); %variables
dydt=[-b*s*i;b*s*i-a*i;a*i];
end

function J=oriSIRD_Jacobian(t,y,params)
a=params(1);b=params(2);
s=y(1);i=y(2);n=y(3); %variables
J=[-b*i -b*s 0;b*i b*s-a 0;0 a 0];
end

function dydt=modSIRD_Vpool(t,y,params)
% modified SIRD model with a variable representing for the pool of virus
c=params(1);d=params(2);k=params(3);lambda=params(4);%n=params(5);b=params(6)  %parameters
S=y(1);I=y(2);V=y(4); %variables
%Hill_term=(V^n)/(V^n+b);
dydt=[-k*S*V;k*S*V-lambda*I;lambda*I;c*I-d*V];
end

function J=modSIRD_Vpool_Jacobian(t,y,params)
% Jacobian matrix for modified SIRD model with a variable representing for the pool of virus
c=params(1);d=params(2);k=params(3);lambda=params(4);  %parameters
S=y(1);V=y(4); %variables
J=[-k*V 0 0 -k*S;k*V -lambda 0 k*S;0 lambda 0 0;0 c 0 -d];
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

function [x_opt,f_opt,x_kept,f_kept]=fitSIRDModel(dydtModelFunc,JacobianModelFunc,...
    model_name,with_virus,data_cases_active,data_cases_total,...
    total_population,params_lb,params_ub,days_prediction,start_date)
%Function for fitting a specific pandemic model (dydtModelFunc,
%JacobianModelFunc) to population data
% Input variables:
% - 1. dydtModelFunc: function handle for the function:
%      dydt=dydtModelFunc(t,y,params) which describes the dynamics of virus
%      spreading.
% - 2. JacobianModelFunc: function handle for the function
%      J=JacobianModelFunc(t,y,params) which describes the Jacobian matrix
%      (i.e. J=dy/d(params)) of the model described in dydtModelFunc
% - 3. model_name: string variable for the name of the model
% - 4. with_virus: label indicating whether virus pool is included in the
%      model. 1 for yes; 0 for no
% - 5. data_infected: numbers of infected peoples on different dates
% - 6. data_not_susceptible: numbers of died/recovered people on different
%      dates
% - 7. total_population: total population in the country studied
% - 8. params_lb: lower bounds for parameters in the model
% - 9. params_ub: upper bounds for parameters in the model
% - 10. start_date: first date with non-zero active case number

data_not_susceptible = data_cases_total - data_cases_active;
days_data_available = 0:length(data_cases_total)-1;

if with_virus
    y_initial=[total_population-data_cases_total(1) data_cases_active(1) data_not_susceptible(1) 0];
else
    y_initial=[total_population-data_cases_total(1) data_cases_active(1) data_not_susceptible(1)];
end
 
DSA_options.start_beta=1e-6;
DSA_options.end_beta=1e8;
DSA_options.expected_fobj=length(days_data_available)*0.00005;
DSA_options.visual_output=0;

[x_opt,f_opt,x_kept,f_kept] = DSA(@calResidues,log10(params_lb(:)),...
    log10(params_ub(:)),DSA_options);

I_samp=zeros(size(x_kept,1),length(days_prediction));
N_samp=zeros(size(x_kept,1),length(days_prediction));
for i=1:size(x_kept,1)
    params=10.^x_kept(i,:);
    opts=odeset('Jacobian',@(t,y)JacobianModelFunc(t,y,params));
    [t,y]=ode23s(@(t,y)dydtModelFunc(t,y,params),days_prediction,y_initial,opts);
    I_samp(i,:)=y(:,2)';
    N_samp(i,:)=y(:,3)';
end
I_lb=min(I_samp);I_ub=max(I_samp);
N_lb=min(N_samp);N_ub=max(N_samp);

params_opt=10.^x_opt;
opts=odeset('Jacobian',@(t,y)JacobianModelFunc(t,y,params_opt));
[t,y]=ode23s(@(t,y)dydtModelFunc(t,y,params_opt),days_prediction,y_initial,opts);

if size(x_kept,1)==0
    I_lb=y(:,2)';
    I_ub=y(:,2)';
    N_lb=y(:,3)';
    N_ub=y(:,3)';
end

figure;
subplot(1,2,1);
fill([t' t(end:-1:1)'],[I_lb I_ub(end:-1:1)],[0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on;plot(t,y(:,2));scatter(days_data_available,data_cases_active,'o');
title(strcat(model_name,':Infected'));
xlabel(sprintf('Days since %s',start_date));ylabel('Count');
subplot(1,2,2);
fill([t' t(end:-1:1)'],[N_lb N_ub(end:-1:1)],[0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on;plot(t,y(:,3));scatter(days_data_available,data_not_susceptible,'o');
title(strcat(model_name,':Recovered/dead'));
xlabel(sprintf('Days since %s',start_date));ylabel('Count');
saveas(gcf,strcat(model_name,'.fig'));

    function f=calResidues(log_params)
        params=10.^log_params;
        tstart=tic;
        
        opts=odeset('Jacobian',@(t,y)JacobianModelFunc(t,y,params),...
            'Events',@(t,y) Integration_Too_Slow(t,y,params,tstart),'RelTol',1e-6,'AbsTol',1e-6);
        [t,y]=ode23s(@(t,y)dydtModelFunc(t,y,params),days_data_available,y_initial,opts);
        if size(y,1)==length(days_data_available)
            f=[(y(:,2)-data_cases_active(:))/max(data_cases_active);(y(:,3)-...
                data_not_susceptible(:))/max(data_not_susceptible)];
        else
            f=10000*ones(2*length(data_cases_active),1);
        end
        if sum(isnan(y(:)))>0
            f=10000*ones(2*length(data_cases_active),1);
        end
        f=real(f);
    end
end
