clear all
close all
suffix = "_Allen";
results_folder = "/Users/bognasmug/Q-NQ-data-analysis/Matlab_scripts/results_2020_11/Allen";
mkdir(results_folder);
options=odeset('RelTol',1e-13);

% experimental 
G0 = 111;
mg_count=33*10^6; %cells in mg
P0=0.01; %biomass [in g] used at the beginning of fresh cells growth
P0starvation = 2*10^7/mg_count;
% regrowth: same parameters as at the initial phase BUT with a lag
regrowth_times = linspace(1,23,24);


% fitted in R using DEOPTIM
Vh = 130;
Kh = 91;
a = 0.007;


%2) other params

% fitted to Allen
dQ= 0.002;
dNQ = 0.005;


% fitted to LeeLeu
%dQ = 0.0001;
%dNQ = 0.003;

reusablility_ratio = 0.5;
% this one will not change regardless of the staving environment;
reusablility_ratio_regrowth = 0.5;

epsilon_NQ = 0.2;
epsilon_Q = 0.25;
epsilon = 0.02;
epsilonAG = 0.0001;

% switching from NQ results ins slightly less Q. This is to represent the
% fact that cells that switch to Q are packed with resources and could have
% potentially duplicated
b= 0.9;

% conditions to span through
freqQ = [0 0.75 1];
types = ['N', 'S', 'Q'];
Sidx = 2; % wjoch one is S


% our lag data
LAGS_H20_Q = [1.5, 1.8 1.7 2.8 3.3 1.6];
LAGS_H20_NQ = [2.2, 4.1 5 7.5 7.63 9.4];
LAGS_YPD_Q = [1.6, 2.1 1.5 2.8 2.8 2.6];
LAGS_YPD_NQ = [2, 2.1 2.38 4 7.2 6.25];


%lags for YPD same as for H20
%LAGS_YPD_Q = LAGS_H20_Q;
%LAGS_YPD_NQ = LAGS_H20_NQ;

figure()
G = linspace(0,60,100);
eps = switch_rate_of_glucose(G);
plot(G, eps*epsilon_NQ, ':r', 'LineWidth', 4)
hold on
plot(G,  eps*epsilon_Q, '-b', 'LineWidth', 4)
legend('descendants of NQ', 'descendants of Q')
ylabel('Switch rate to Q (\epsilon^{QQ}, \epsilon^{NQ} )')
xlabel('Glucose concentration[mM]')


%% very short starvation (4 days) results in ~ 75% Q cells when the culture is started from NQ descendants
scenario = 'culture_differentiation_within_4days';
time_growth = 4*24;
y0SNQ = [G0; 0; 0; P0;0;0;0];
[tNQ,YsNQ] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, time_growth], y0SNQ, options);    
Qcells = YsNQ(:,2);
NQcells = YsNQ(:,3)+YsNQ(:,4) +YsNQ(:,5);
figure()
hold on
plot(tNQ, YsNQ(:,2), '-b')
plot(tNQ, YsNQ(:,3)+YsNQ(:,4) +YsNQ(:,5), '-r')
plot(tNQ, YsNQ(:,2)+YsNQ(:,3)+YsNQ(:,4) +YsNQ(:,5), '-k')
legend('Q', 'NQ', 'total')
xlabel('time [h]')
ylabel('biomass: NQ descendants')
title(scenario)

varNames = {'Growth_Time', 'Q_biomass', 'NQ_biomass'};
dataT=table(tNQ, Qcells, NQcells, 'VariableNames',varNames);
writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))


%% Starvation in water: various weeks
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20';

% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants

for week = 1:6
   Qlag = LAGS_H20_Q(week);
   NQlag = LAGS_H20_NQ(week);
   Tstarvation = week*7*24;
   for i = 1:length(freqQ)
        f=freqQ(i); 
        type = types(i);
        for j = 1:length(regrowth_times)
            % start starvation with some Q and some NQ
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation H20 (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))

%% Starvation in water: various weeks: assuming Q lags are the same from week to week
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20_assuming_same_Q_lags';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants

for week = 1:6
   Qlag = mean(LAGS_H20_Q);
   NQlag = LAGS_H20_NQ(week); 
   Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
    f=freqQ(i); 
    type = types(i);
    for j = 1:length(regrowth_times)
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));

            %{
            figure()
            hold on
            plot(t, Ys(:,2), '-b')
            plot(t, Ys(:,3)+Ys(:,4) +Ys(:,5), '-r')
            plot(t, Ys(:,2)+Ys(:,3)+Ys(:,4) +Ys(:,5), '-k')
            legend('Q', 'NQ', 'total')
            xlabel('time [h]')
            ylabel('biomass: NQ descendants')
            %}
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))
%% Starvation in water: various weeks: 
% assuming no death rates
% lags varying with week times
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20_assuming_same_no_death_rate_and_constant_Qlag';

% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
dQ0 = 0;
dNQ0 = 0;
for week = 1:6 
   Qlag = mean(LAGS_H20_Q);
   %Qlag = LAGS_H20_Q(week);
   NQlag =LAGS_H20_NQ(week);
   Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
    f=freqQ(i); 
    type = types(i);
    for j = 1:length(regrowth_times)
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ0, dNQ0, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ0, dNQ0, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))
%% Starvation in water: various weeks: 
% assuming Q lags equal NQ lags
% and dQ equals to dNQ
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20_assuming_Qlags_equal_NQlags';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
for week = 1:6
   NQlag = LAGS_H20_NQ(week);
   Qlag = NQlag; 
   Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
    f=freqQ(i); 
    type = types(i);
    for j = 1:length(regrowth_times)
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dNQ, dNQ, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dNQ, dNQ, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))

%% Starvation in water: various weeks: 
% assuming no death rates
% lags constant Q < NQ
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20_assuming_same_no_death_rate_constant_Qlag_smaller_than_NQ_lag';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
dQ0 = 0;
dNQ0 = 0;
for week = 1:6 
   Qlag = mean(LAGS_H20_Q);
   NQlag = mean(LAGS_H20_NQ); %LAGS_H20_NQ(week); 
   Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
    f=freqQ(i); 
    type = types(i);
    for j = 1:length(regrowth_times)
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ0, dNQ0, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ0, dNQ0, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/',scenario, suffix, '.txt'))

%% Starvation in water: various weeks: 
% assuming zero lags
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_H20_assuming_zero_lag';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
for week = 1:6
   Qlag = 0;
   NQlag = 0; 
   Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
    f=freqQ(i); 
    type = types(i);
    for j = 1:length(regrowth_times)
            Q0=f*P0starvation;
            NQ_NQ0=(1-f)*P0starvation;
            y0S = [0;Q0;0;NQ_NQ0;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);          
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
    end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
    subplot(6,1,week)
    title(strcat('Long starvation (',num2str(week), ' week)'))
    hold on
    for i = 1:length(freqQ)
    plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
    end
    xlabel('regrowth time [h]')
    ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
    legend()
    title(scenario)
    ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))

%% Starvation in YPD: multiple weeks
% assume in regrowth there is no differentiation to NQ and Q
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type'};
dataT=array2table(zeros(0,4), 'VariableNames',varNames);
scenario = 'long_starvation_YPD';

% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
reusablility_ratio = 0.5; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
freqQ = [0 0.75 1];
Sidx = 2;
for week = 1:6
    Qlag = LAGS_YPD_Q(week);
    NQlag = LAGS_YPD_NQ(week);
    Tstarvation = week*7*24;
  for i = 1:length(freqQ)
   % start starvation with some Q and some NQ
   f=freqQ(i); 
   type = types(i);
     for j = 1:length(regrowth_times) 
            y0S = [0; f*P0starvation; 0; (1-f)*P0starvation;0;0;0];
            [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    
            % proportion of Q and NQ after starvation
            final_prop_long_starvation(i,j, week)=Ys(end,2)/(Ys(end,2) + Ys(end,3)  + Ys(end,4)  + Ys(end,5));
            % now the regrowth
            y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);  
            % final amount of biomass
            final_biomass_long_starvation(i,j, week) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
end
end

S=final_biomass_long_starvation(Sidx,:,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_long_starvation = final_biomass_long_starvation./S;

figure()
for week=1:6
subplot(6,1,week)
  title(strcat('Long starvation in YPD (',num2str(week), ' week'))
hold on
for i = 1:length(freqQ)
plot(regrowth_times, relative_fitness_to_S_long_starvation(i,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
end
xlabel('regrowth time [h]')
ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
legend()
    title(scenario)
ylim([0 1.5])
end

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))

%% Starvation: very short in YPD i.e. 4 days
% show that not switching to Q may be actually beneficial

% assume in regrowth there is no differentiation to NQ and Q
varNames = {'Regrowth_Time', 'finalBiomass',  'type'};
dataT=array2table(zeros(0,3), 'VariableNames',varNames);
scenario = 'short_starvation_scenario_YPD';
Qlag = 2;
NQlag = 1;

final_biomass_short_starvation = [];
reusablility_ratio = 0.5; % that means G=0 -> switching rate =0 -> not important if our NQ are Qdescendants or NQdescendants
Sidx = 2;
  for i = 1:length(freqQ)
   % start starvation with some Q and some NQ
   f=freqQ(i); 
   type = types(i);
     for j = 1:length(regrowth_times) 
         % now the regrowth
            y0R = [G0;0; f*P0starvation; (1-f)*P0starvation;0;0;0];
            [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);  
         % final amount of biomass
            final_biomass_short_starvation(i,j) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
     end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_short_starvation(i,:)';
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
   end

S=final_biomass_short_starvation(Sidx,:);
% (N1end/N10)/(N2end/N20) = N1end/N2end as here they start from same N0
relative_fitness_to_S_short_starvation = final_biomass_short_starvation./S;

figure()
  title(strcat('short starvation in YPD'))
hold on
for i = 1:length(freqQ)
plot(regrowth_times, relative_fitness_to_S_short_starvation(i,:), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
end
xlabel('regrowth time [h]')
ylabel({'relative fitness of the population', 'to the one where f_Q=75%'})
legend()
    title(scenario)
ylim([0 1.5])

writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))


%% Compare NQ in H20 and YPD
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type', 'environment'};
dataT=array2table(zeros(0,5), 'VariableNames',varNames);
scenario = 'long_starvation_NQ_YPD_VS_H20';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
for env = [1 2]
    if env == 1
        environment = "H20";
        reusablility_ratio=0;
        LAGSQ = LAGS_H20_Q;
        LAGSNQ=LAGS_H20_NQ;
    else
        environment = "YPD";
        reusablility_ratio=0.5;
        LAGSQ = LAGS_YPD_Q;
        LAGSNQ=LAGS_YPD_NQ;

    end
for week = 1:6
    Qlag = LAGSQ(week);
    NQlag = LAGSNQ(week);
    Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
        f=freqQ(i); 
        type = types(i);
        for j = 1:length(regrowth_times)
                y0S = [0; f*P0starvation; 0; (1-f)*P0starvation;0;0;0];
                [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    

                % now the regrowth
                y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
                [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);  
         
                %[t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, 0, 0, 0,0, 0, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);          

                % final amount of biomass
                final_biomass_long_starvation(i,j, week, env) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week, env)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataTnew.environment = repmat(environment, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
end
end
end

NQYPDtoNQH20 = final_biomass_long_starvation(1,:,:,2)./final_biomass_long_starvation(1,:,:,1);
figure()
for week=1:6
subplot(6,1,week)
title('Long starvation NQ in H20 vs YPD comparison')
hold on
plot(regrowth_times, NQYPDtoNQH20(1,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
xlabel('regrowth time [h]')
ylabel({'relative fitness of NQ in YPD to H20'})
title(scenario)
end


writetable(dataT,strcat(results_folder, '/',scenario, suffix, '.txt'))


%% additional checks
% check that what makes NQ do better in YPD is the reincarnation
% assume lags in both scenarios are the same (zero) 
varNames = {'Regrowth_Time', 'finalBiomass', 'week', 'type', 'environment'};
dataT=array2table(zeros(0,5), 'VariableNames',varNames);
scenario = 'long_starvation_NQ_YPD_VS_H20_assuming_no_Lags';
% Experiment 1a
final_prop_long_starvation = [];
final_biomass_long_starvation = [];
for env = [1 2]
    if env == 1
        environment = "H20";
        reusablility_ratio=0;

    else
        environment = "YPD";
        reusablility_ratio=0.5;

    end
for week = 1:6
    Qlag = 0;
    NQlag = 0;
    Tstarvation = week*7*24;
    for i = 1:length(freqQ)
    % start starvation with some Q and some NQ
        f=freqQ(i); 
        type = types(i);
        for j = 1:length(regrowth_times)
                y0S = [0; f*P0starvation; 0; (1-f)*P0starvation;0;0;0];
                [t,Ys] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, 0, 0, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b),[0, Tstarvation], y0S, options);    

                % now the regrowth
                y0R = [G0; 0; Ys(end,2);Ys(end,3)+Ys(end,4)+Ys(end,5);0;0;0];
                %[t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, 0, 0, 0,0, 0, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options); 
                [t,Y] = ode15s(@(t,y) simulate_model(t,y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio_regrowth, dQ, dNQ, b),[0, regrowth_times(j)], y0R, options);  
            

                % final amount of biomass
                final_biomass_long_starvation(i,j, week, env) = Y(end,2) + Y(end,3)  + Y(end,4)  + Y(end,5);
        end
        dataTnew=table(regrowth_times', 'VariableNames',{'Regrowth_Time'});
        dataTnew.finalBiomass = final_biomass_long_starvation(i,:, week, env)';
        dataTnew.week = repmat(week, length(regrowth_times),1);
        dataTnew.type = repmat(type, length(regrowth_times),1);
        dataTnew.environment = repmat(environment, length(regrowth_times),1);
        dataT = [dataT;dataTnew];
end
end
end

NQYPDtoNQH20 = final_biomass_long_starvation(1,:,:,2)./final_biomass_long_starvation(1,:,:,1);

figure()
for week=1:6
subplot(6,1,week)
title('Long starvation NQ in H20 vs YPD, assuming same lags')
hold on
plot(regrowth_times, NQYPDtoNQH20(1,:, week), 'DisplayName', strcat('f_Q = ',num2str(freqQ(i))))
xlabel('regrowth time [h]')
ylabel({'relative fitness of NQ in YPD to H20'})
title(scenario)
end


writetable(dataT,strcat(results_folder, '/', scenario, suffix, '.txt'))
