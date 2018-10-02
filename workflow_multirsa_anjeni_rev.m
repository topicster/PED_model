%% 1. REGIONAL SENSITIVITY ANALYSIS OF THE PED MODEL IN THE ANJENI WATERSHED
disp('REGIONAL SENSITIVITY ANALYSIS OF THE PED MODEL IN THE ANJENI WATERSHED')

% This script provides an application example of Regional Sensitivity
% Analysis (RSA) (see help RSA_indices_thres.m for more details about this
% Uncertainty Analysis method and references).
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff-erosion PED model developed
% for the Ethiopian highlands by Collick et al. (2009), Steenhuis et al.
% (2009), Tesemma et al. (2010), and Tilahun et al. (2013a, 2013b, 2015).
%
% The model is applied to the Anjeni catchment in Ethiopia, Africa
% (see header of file Anjeni_catchment.txt for more details).
% The inputs subject to SA are the 9 model parameters, and the scalar 
% output for SA is (one or multiple) performance metric.
% Alternatively, the inputs subject to SA are the 13 coupled
% hydrological-erosion model parameters.
%
% A different setup runs the model hymod with the same data.
%
% The resulting FLUXES from the model best prediction (or other asigned by
% the user) are input in an erosion model by Tilahun et al. (2013a, 2013b,
% 2015) that computes sediment concentration.
% The inputs for the erosion model are the 4 model parameters, a time
% series of runoff-producing area with active rill formation, and the
% FLUXES consisting of runoff in the 2 areas, baseflow, and interflow.
%
% This script has been edited by Boris Ochoa-Tocachi in 2017 based on a
% script prepared by Francesca Pianosi and Fanny Sarrazin (2014).
% mail to: boris.ochoa13@imperial.ac.uk

disp(' ')
%% 2. INITIALISE THE DIRECTORY
disp('INITIALISE THE DIRECTORY')

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

disp(' ')
%% 3. INITIALISE THE DATA
disp('INITIALISE THE DATA')

% Load data:
data = load('Anjeni_catchment.txt','-ascii');
areakm2 = 1.13; % km^2

% [T,~] = size(data);
% rain = data(:,1);
% evap = data(:,2);
% flow = data(:,3);
% sedm = data(:,4);
% hvec = data(:,5);

% Define the data range in years
yearini = 1990;
yearend = 1993;

% Approximate (no leap years)
% ini = 365*(yearini-1984)+1;
% fin = 365*(yearend-1984+1);

% With lap years:
% Extract data from 1990 to 1993
ini = 2193; % 01 Jan 1990
fin = 3653; % 31 Dec 1993

T = fin-ini+1;
rain = data(ini:fin,1);
evap = data(ini:fin,2);
flow = data(ini:fin,3);
sedm = data(ini:fin,4);

% Define output function:
% myfun = 'ethiopian_sim' ;
% myfun = 'ethiopian_mod' ;
myfun = 'ethiopian_same' ;
warmup = 90 ; % warmup period to be discarded

% Plot data:
figure
subplot(311)
plot(rain)
ylabel('Rain (mm day^{-1})')
set(gca,'XLim',[0 T+1])
subplot(312)
plot(evap)
ylabel('PET (mm day^{-1})')
set(gca,'XLim',[0 T+1])
subplot(313)
plot(flow)
ylabel('Flow (mm day^{-1})')
xlabel('time (days)')
set(gca,'XLim',[0 T+1])
title('mm day^{-1}')

disp(' ')
%% 4. HYDROLOGICAL EXPERIMENT SETUP
disp('HYDROLOGICAL EXPERIMENT SETUP')

M  = 9 ; % number of uncertain parameters: 9
% 1. A1  = Percentage of saturated area [%]
% 2. A2  = Percentage of degraded area [%]
% 3. A3  = Percentage of infiltration area [%]
%          A1+A2+A3 = 100 (%) (see line 82)
% 4. Sm1 = Maximum water storage A1 [mm]
% 5. Sm2 = Maximum water storage A2 < A1 [mm]
% 6. Sm3 = Field capacity of A3, generally > A2 > A1 [mm]
% 7. BSm = Maximum baseflow storage of the lineal reservoir [mm]
% 8. tm  = half-life time of the baseflow linear reservoir [day]
% 9. tau = duration of the interflow zero-order reservoir [day]
DistrFun  = 'unif'  ; % Parameter distribution
% Parameter ranges:
%        A1  A2  A3 Sm1 Sm2 Sm3 BSm t(1/2) tau*
%      [  2  10  47 200  10  65 100    75   10]; % Guzman et al., 2017.
%      [  2  14  50 200  10 100 100    70   10]; % Tilahun et al., 2013.
%      [  0  20  60   - 150 250  70    70   20]; % Engda et al., 2012; Legesse et al. 2009
xmin = [  0   5  20  50   5  50  50    20    1] ; % minimum values
xmax = [ 20  40  80 250 150 300 200   100   50] ; % maximum values

% Parameter names
X_Labels = {'A_1','A_2','A_3','S_{max_1}','S_{max_2}','S_{max_3}','BS_{max}','t_{1/2}','\tau^{*}'} ;

disp(' ')
% SAMPLING INPUT SPACE
disp('SAMPLING INPUT SPACE')

% Sample the input space using the 'AAT_sampling' function
% (see help AAT_sampling.m for more details)
SampStrategy = 'lhs' ; % sampling strategy (other options is 'rsu')
% Latin Hypercube Sampling
N = 10000 ; % sample size
DistrPar = cell(M,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:M; DistrPar{i} = [xmin(i) xmax(i) ] ; end % this is to resort to
% the format required by AAT_sampling function.
% Perform sampling:
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);
% Normalise the areas
% counter = 0; % Just for information
for i = 1:N
    if sum(X(i,1:3)) > 100
        % counter = counter + 1;
        X(i,1:3) = 100*X(i,1:3)/sum(X(i,1:3));
    end
end

% Plot results of the random sampling
figure; parcoor(X,X_Labels); % Parallel coordinate plots

disp(' ')
%% 5. MODEL EVALUATION
disp('MODEL EVALUATION')

% Run the model and compute model output at sampled parameter sets:
Qsim = model_evaluation(myfun,X,rain,evap);

% Plot simulated time series:
figure
subplot(6,1,1)
plot(rain,'k')
ylabel('Rain (mm day^{-1})')
set(gca,'XLim',[0 T+1])
subplot(6,1,2)
plot(evap,'k')
ylabel('PET (mm day^{-1})')
set(gca,'XLim',[0 T+1])
% subplot(6,1,3)
% plot(PETsim,'k')
% ylabel('PET (mm day^{-1})')
% set(gca,'XLim',[0 T+1])
subplot(6,1,3:6)
a = plot(Qsim','Color',[0.33 0.33 0.33]);
% add flow observations as black circles
hold on
c = plot(flow,'o','Color',[0 0 0],'MarkerSize',2);
% Customise the plot
set(gca,'XLim',[0 T+1],'YLim',[-2.5 50])
xlabel('time (days)'); ylabel('Flow (mm day^{-1})')
box on
grid on
legend([a(1) c],'Simulated flow','Observed')
clear a c

disp(' ')
%% 6. MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS
disp('MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS')

% Statistics
hydro_nse = NSE(Qsim(:,warmup+1:end),flow(warmup+1:end)');
hydro_pbias = PBIAS(Qsim(:,warmup+1:end),flow(warmup+1:end)');
hydro_kge = KGE(Qsim(:,warmup+1:end),flow(warmup+1:end)');
% Transformed objectives for minimisation
hydro_NSEmod = 1 - (hydro_nse); % 1-NSE to minimise
hydro_PBIASmod = abs(hydro_pbias); % abs(PBIAS) to minimise
hydro_KGEmod = 1 - (hydro_kge); % 1-KGE to minimise

% Visualize input/output samples (this may help finding a reasonable value
% for the output threshold):
figure; scatter_plots(X,hydro_nse,[],'NSE',X_Labels);
figure; scatter_plots(X,hydro_pbias,[],'PBIAS',X_Labels);
figure; scatter_plots(X,hydro_kge,[],'KGE',X_Labels);

% Set user-defined output thresholds
% nse_thres = 0.7; % behavioural threshold for NSE
% pbias_thres = 5; % behavioural threshold for PBIAS
% kge_thres = 0.8; % % behavioural threshold for KGE, or same as NSE

% OR: Define the thresholds from the best values of each Obj
hydro_nse_sortindx = []; hydro_pbias_sortindx = []; hydro_kge_sortindx = [];
[hydro_nse_sortindx(:,1),hydro_nse_sortindx(:,2)] = sort(hydro_NSEmod);
[hydro_pbias_sortindx(:,1),hydro_pbias_sortindx(:,2)] = sort(hydro_PBIASmod);
[hydro_kge_sortindx(:,1),hydro_kge_sortindx(:,2)] = sort(hydro_KGEmod);
% Extract thresholds
hydro_nse_thres = 1-hydro_nse_sortindx(ceil(N*0.05));
hydro_pbias_thres = abs(hydro_pbias_sortindx(ceil(N*0.05)));
hydro_kge_thres = 1-hydro_kge_sortindx(ceil(N*0.05));

% Best fit using criteria:
% NSE
hydro_nse_sortset = X(hydro_nse_sortindx(:,2),:);
hydro_nse_bestObj = [hydro_nse(hydro_nse_sortindx(1,2),:) hydro_pbias(hydro_nse_sortindx(1,2),:) hydro_kge(hydro_nse_sortindx(1,2),:)];
% PBIAS
hydro_pbias_sortset = X(hydro_pbias_sortindx(:,2),:);
hydro_pbias_bestObj = [hydro_nse(hydro_pbias_sortindx(1,2),:) hydro_pbias(hydro_pbias_sortindx(1,2),:) hydro_kge(hydro_pbias_sortindx(1,2),:)];
% KGE
hydro_kge_sortset = X(hydro_kge_sortindx(:,2),:);
hydro_kge_bestObj = [hydro_nse(hydro_kge_sortindx(1,2),:) hydro_pbias(hydro_kge_sortindx(1,2),:) hydro_kge(hydro_kge_sortindx(1,2),:)];
% MULTI
% Normalise the objectives
hydro_NSEnorm   = (hydro_NSEmod-min(hydro_NSEmod))/(max(hydro_NSEmod)-min(hydro_NSEmod));
hydro_PBIASnorm = (hydro_PBIASmod-min(hydro_PBIASmod))/(max(hydro_PBIASmod)-min(hydro_PBIASmod));
hydro_KGEnorm = (hydro_KGEmod-min(hydro_KGEmod))/(max(hydro_KGEmod)-min(hydro_KGEmod));
hydro_multi = sqrt(hydro_NSEnorm.^2 + hydro_PBIASnorm.^2 + hydro_KGEnorm.^2);
% Best fit using Multi criteria
hydro_multi_sortindx = [];
[hydro_multi_sortindx(:,1),hydro_multi_sortindx(:,2)] = sort(hydro_multi);
hydro_multi_sortset = X(hydro_multi_sortindx(:,2),:);
hydro_multi_bestObj = [hydro_nse(hydro_multi_sortindx(1,2),:) hydro_pbias(hydro_multi_sortindx(1,2),:) hydro_kge(hydro_multi_sortindx(1,2),:)];

% Simulations that comply with Moriasi et al. (2007) recommendations
hydro_M2007 = and(hydro_NSEmod<0.5,hydro_PBIASmod<25); % location
hydro_M2007n = 100*sum(hydro_M2007)/N; % percentage

disp(' ')
%% 7. RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS
disp('RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS')

% RSA (find behavioural parameterizations):
hydro_threshold = [1-hydro_nse_thres hydro_pbias_thres 1-hydro_kge_thres];
[hydro_multi_mvd_rsa,hydro_multi_idxb_rsa] = RSA_indices_thres(X,[hydro_NSEmod hydro_PBIASmod hydro_KGEmod],hydro_threshold);
hydro_multi_Xbeh = X(hydro_multi_idxb_rsa,:);

hydro_top5n = 100*sum(hydro_multi_idxb_rsa)/N; % percentage of top 5% behavioural

% Highlight the behavioural parameterizations in the scatter plots:
figure; scatter_plots(X,hydro_nse,[],'NSE',X_Labels,hydro_multi_idxb_rsa);
for i=1:M
    subplot(ceil(M/5),5,i)
    plot(hydro_nse_sortset(1,i),hydro_nse_bestObj(1),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(hydro_pbias_sortset(1,i),hydro_pbias_bestObj(1),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(hydro_kge_sortset(1,i),hydro_kge_bestObj(1),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(hydro_multi_sortset(1,i),hydro_multi_bestObj(1),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set')

figure; scatter_plots(X,hydro_pbias,[],'PBIAS',X_Labels,hydro_multi_idxb_rsa);
for i=1:M
    subplot(ceil(M/5),5,i)
    plot(hydro_nse_sortset(1,i),hydro_nse_bestObj(2),'vk','MarkerSize',10,'MarkerFaceColor',0.6*[1 1 1])
    plot(hydro_pbias_sortset(1,i),hydro_pbias_bestObj(2),'^k','MarkerSize',10,'MarkerFaceColor',0.7*[1 1 1])
    plot(hydro_kge_sortset(1,i),hydro_kge_bestObj(2),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(hydro_multi_sortset(1,i),hydro_multi_bestObj(2),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set')

figure; scatter_plots(X,hydro_kge,[],'KGE',X_Labels,hydro_multi_idxb_rsa);
for i=1:M
    subplot(ceil(M/5),5,i)
    plot(hydro_nse_sortset(1,i),hydro_nse_bestObj(3),'vk','MarkerSize',10,'MarkerFaceColor',0.6*[1 1 1])
    plot(hydro_pbias_sortset(1,i),hydro_pbias_bestObj(3),'^k','MarkerSize',10,'MarkerFaceColor',0.7*[1 1 1])
    plot(hydro_kge_sortset(1,i),hydro_kge_bestObj(3),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(hydro_multi_sortset(1,i),hydro_multi_bestObj(3),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')

% Plot parameter CDFs:
RSA_plot_thres(X,hydro_multi_idxb_rsa,[],X_Labels,{'behav','non-behav'}); % add legend

disp(' ')
%% 8. SENSITIVITY INDICES (MVD) USING NSE, PBIAS, AND KGE
disp('SENSITIVITY INDICES (MVD) USING NSE, PBIAS, AND KGE')

% mvd: maximum vertical distance between parameters CDFs

% Assess robustness by bootstrapping:
Nboot = 100;
[hydro_multi_mvd,hydro_multi_idxb,hydro_multi_mvd_lb,hydro_multi_mvd_ub] = RSA_indices_thres(X,[hydro_NSEmod hydro_PBIASmod hydro_KGEmod],hydro_threshold,1,Nboot);
[hydro_multi_spread,~,hydro_multi_spread_lb,hydro_multi_spread_ub] = RSA_indices_thres(X,[hydro_NSEmod hydro_PBIASmod hydro_KGEmod],hydro_threshold,2,Nboot);
% Plot results:
figure; boxplot1(hydro_multi_mvd,X_Labels,'mvd',hydro_multi_mvd_lb,hydro_multi_mvd_ub)
figure; boxplot1(hydro_multi_spread,X_Labels,'spread',hydro_multi_spread_lb,hydro_multi_spread_ub)
set(gca,'YLimMode','auto')

disp(' ')
%% 9. SET UP OF THE EROSION MODEL
disp('SET UP OF THE EROSION MODEL')

% Non-negative sediment values are discarded
% sedm(~(sedm>0)) = nan;
% sedm = data(ini:fin,4);

% Fraction of the runoff-producing area with active rill formation.
% H = 1.00; % middle of June to middle of July
% H = 0.50; % middle of July to first week of August
% H = 0.25; % first to second week of August
% H = 0.00; % remainder season

% Discrete H for a year of 365 days (Tilahun et al., 2013a, 2013b)
H_disc = zeros(365,1);  % All year
H_disc(167:196) = 1.00; % 16 Jun - 15 Jul
H_disc(197:212) = 0.50; % 16 Jul - 31 Jul
H_disc(213:227) = 0.25; % 01 Aug - 15 Aug

% Continuous H for a year of 365 days (based on Guzman et al., 2017)
H_pnt = NaN(365,1);
% H_pnt(1)   = 0.00; % 01 Jan
H_pnt(66)  = 0.00; % 01 Mar
% H_pnt(91)  = 0.00;   % 01 Apr
H_pnt(128) = 0.25; % 08 Apr
% H_pnt(136) = 0.25;   % 16 May
H_pnt(152) = 0.60;   % 01 Jun
H_pnt(182) = 1.00;   % 01 Jul
H_pnt(212) = 0.50;   % 31 Jul
H_pnt(227) = 0.25;   % 15 Aug
H_pnt(273) = 0.00;   % 31 Sep
% H_pnt(365) = 0.00; % 31 Dec
% Create a curve using cubic spline interpolation
valid_pnt = find(~isnan(H_pnt));
H_csi = csape(valid_pnt,H_pnt(valid_pnt),'periodic');
H_fval = fnval(H_csi,66:273);
H_cont = [zeros(65,1);H_fval';zeros(365-274+1,1)];
H_cont(H_cont<0) = 0; % Correct any negative values
H_cont(H_cont>1) = 1; % Correct values above one

% Repeat H for the analysed time interval
% H = repmat(H_cont,(yearend-yearini+1),1);
% Create H from 1984 to 1993 with leap years
%    1984 - 1987          1988 - 1991          1992 - 1993
hvec = [0;repmat(H_cont,4,1);0;repmat(H_cont,4,1);0;repmat(H_cont,2,1)];
hvec2 = [0;repmat(H_disc,4,1);0;repmat(H_disc,4,1);0;repmat(H_disc,2,1)];
% Extract data from 1990 to 1993
H = hvec(ini:fin);
H2 = hvec2(ini:fin);

% Plot H
figure
% H parameter
a = plot(H_disc,'-','Color',[0.66 0.66 0.66],'Linewidth',2);
hold on
plot(H_pnt,'ok')
b = plot(H_cont,'-.','Color',[0.33 0.33 0.33],'Linewidth',1);
% Cumulative curves
H_cum1 = cumsum(H_disc,'omitnan');
plot(H_cum1/max(H_cum1),'-','Color',[0.66 0.66 0.66],'Linewidth',2)
% H_cum2 = cumsum(H_pnt,'omitnan');
% H_cum2(isnan(H_pnt)) = nan;
% plot(H_cum2/max(H_cum2),'ok')
H_cum3 = cumsum(H_cont,'omitnan');
plot(H_cum3/max(H_cum3),'-.','Color',[0.33 0.33 0.33],'Linewidth',1)
% Customise the plot
title('fraction of area with active rill formation')
ylabel('Fraction H (between a_t and a_s)')
set(gca,'XLim',[1 365],'YLim',[0 1.25],'XTick',(1:21:365),'YTick',(0:0.25:1.25))
datetick('x','dd-mmm','keeplimits','keepticks')
legend([a b],'Tilahun et al., 2013; 2015','based on Guzman et al., 2017')
clear a b

disp(' ')
%% 10. EROSION EXPERIMENT SETUP
disp('EROSION EXPERIMENT SETUP')

E  = 4; % number of erosion parameters: 4
% 1. a1t  = Sediment transport limit for A1 [(g/L)/(mm/d)^-0.4]
% 2. a2t  = Sediment transport limit for A2 [(g/L)/(mm/d)^-0.4]
% 3. a1s  = Sediment source limit for A1 [(g/L)/(mm/d)^-0.4]
% 4. a2s  = Sediment source limit for A2 [(g/L)/(mm/d)^-0.4]
% DistrFun  = 'unif'  ; % Parameter distribution
% Parameter ranges:
%       at1 at2 as1 as2
%      [  6   5 5.5 4.5]; % Guzman et al. (2017)
%      [0.2 3.4 0.2 3.4]; % Tilahun et al. (2013a)
%      [  4   4   3   3]; % Tilahun et al. (2013b)
%       (at1-as1) (at2-as2)   as1   as2
ymin = [       0         0      0     0] ; % minimum values
ymax = [       5         5     10    10] ; % maximum values

% Parameter names
Y_Labels = {'a_{t_1}','a_{t_2}','a_{s_1}','a_{s_2}'};
XY_Labels = [X_Labels,Y_Labels];

disp(' ')
% SAMPLING INPUT SPACE
disp('SAMPLING INPUT SPACE')

% Sample the input space using the 'AAT_sampling' function
% (see help AAT_sampling.m for more details)
% SampStrategy = 'lhs' ; % sampling strategy (other options is 'rsu')
% Latin Hypercube Sampling
Nbeh = size(hydro_multi_Xbeh,1);
N2 = ceil(10000/Nbeh); % sample size
% N2 = 250;
DistrPar = cell(E,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:E; DistrPar{i} = [ymin(i) ymax(i)]; end % this is to resort to
% the format required by AAT_sampling function.
% Perform sampling:
Y = AAT_sampling(SampStrategy,E,DistrFun,DistrPar,N2);
% Adjust the parameter values: at = (at-as)+as
Y(:,1) = sum(Y(:,[1 3]),2);
Y(:,2) = sum(Y(:,[2 4]),2);

% Combine behavioural hydro parameters with erosion random parameters
N3 = Nbeh*N2;
XY = nan(N3,M+E);
XY(:,M+1:M+E) = repmat(Y,Nbeh,1);
for i = 1:Nbeh
    XY(N2*(i-1)+1:N2*i,1:M) = repmat(hydro_multi_Xbeh(i,:),N2,1);
end
Y = XY(:,M+1:M+E);

% Plot results of the random sampling
figure; parcoor(XY,XY_Labels); % Parallel coordinate plots

disp(' ')
%% 11. EROSION MODEL EVALUATION
disp('EROSION MODEL EVALUATION')

% Erosion model:
myfun2 = 'ethiopian_ped';

% Best performing set
% input_X = multi_sortset(1,:);
% input_X = [  2  10  47 200  10  65  100    75   10]; % Guzman et al., 2017

% Run the model and compute model output at sampled parameter sets:
% Qero = model_evaluation(myfun,input_X,rain,evap);
[Csim,QCsim,CFLUXES] = model_evaluation(myfun2,XY,H,rain,evap);

% Only consider those values when sediment observations are present
Csim(:,~(sedm>0)) = 0;

% Extract fluxes and states:
Qero = nan(size(Csim));
for i = 1:N
    Qero(i,:)     = QCsim{i}(:,1)';
end

% Plot simulated time series:
figure
subplot(6,1,1)
plot(rain,'k')
ylabel('Rain (mm day^{-1})')
set(gca,'XLim',[0 T+1])
subplot(6,1,2)
plot(evap,'k')
ylabel('PET (mm day^{-1})')
set(gca,'XLim',[0 T+1])
subplot(6,1,3)
a = plot(Qero','Color',[0.33 0.33 0.33]);
% add flow observations as black circles
hold on
c = plot(flow,'-','Color',[0 0 0],'MarkerSize',2);
% Customise the plot
set(gca,'XLim',[0 T+1],'YLim',[-2.5 50])
xlabel('time (days)'); ylabel('Flow (mm day^{-1})')
box on
grid on
legend([a(1) c],'Simulated flow','Observed')
clear a c
subplot(6,1,4:6)
a = plot(Csim','Color',[0.33 0.33 0.33]);
% add sediment observations as black circles
hold on
c = plot(sedm,'-o','Color',[0 0 0],'MarkerSize',2);
% Customise the plot
set(gca,'XLim',[0 T+1],'YLim',[-2.5 100])
xlabel('time (days)'); ylabel('Concentration(mm day^{-1})')
box on
grid on
legend([a(1) c],'Simulated concentration','Observed')
clear a c

disp(' ')
%% 12. MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS
disp('MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS')

% Statistics
erosion_nse = NSE(Csim(:,warmup+1:end),sedm(warmup+1:end)');
erosion_pbias = PBIAS(Csim(:,warmup+1:end),sedm(warmup+1:end)');
erosion_kge = KGE(Csim(:,warmup+1:end),sedm(warmup+1:end)');

% Transformed objectives for minimisation
erosion_NSEmod = 1 - (erosion_nse); % 1-NSE to minimise
erosion_PBIASmod = abs(erosion_pbias); % abs(PBIAS) to minimise
erosion_KGEmod = 1 - (erosion_kge); % 1-KGE to minimise

% Visualize input/output samples (this may help finding a reasonable value
% for the output threshold):
figure; scatter_plots(Y,erosion_nse,[],'NSE',Y_Labels);
figure; scatter_plots(Y,erosion_pbias,[],'PBIAS',Y_Labels);
figure; scatter_plots(Y,erosion_kge,[],'KGE',Y_Labels);

% Set user-defined output thresholds
% erosion_nse_thres = 0.7; % behavioural threshold for NSE
% erosion_pbias_thres = 5; % behavioural threshold for PBIAS
% erosion_kge_thres = 0.8; % % behavioural threshold for KGE, or same as NSE

% OR: Define the thresholds from the best values of each Obj
erosion_nse_sortindx = []; erosion_pbias_sortindx = []; erosion_kge_sortindx = [];
[erosion_nse_sortindx(:,1),erosion_nse_sortindx(:,2)] = sort(erosion_NSEmod);
[erosion_pbias_sortindx(:,1),erosion_pbias_sortindx(:,2)] = sort(erosion_PBIASmod);
[erosion_kge_sortindx(:,1),erosion_kge_sortindx(:,2)] = sort(erosion_KGEmod);
% Extract thresholds
erosion_nse_thres = 1-erosion_nse_sortindx(ceil(N3*0.05));
erosion_pbias_thres = abs(erosion_pbias_sortindx(ceil(N3*0.05)));
erosion_kge_thres = 1-erosion_kge_sortindx(ceil(N3*0.05));

% Best fit using criteria:
% NSE
erosion_nse_sortset = XY(erosion_nse_sortindx(:,2),:);
erosion_nse_bestObj = [erosion_nse(erosion_nse_sortindx(1,2),:) erosion_pbias(erosion_nse_sortindx(1,2),:) erosion_kge(erosion_nse_sortindx(1,2),:)];
% PBIAS
erosion_pbias_sortset = XY(erosion_pbias_sortindx(:,2),:);
erosion_pbias_bestObj = [erosion_nse(erosion_pbias_sortindx(1,2),:) erosion_pbias(erosion_pbias_sortindx(1,2),:) erosion_kge(erosion_pbias_sortindx(1,2),:)];
% KGE
erosion_kge_sortset = XY(erosion_kge_sortindx(:,2),:);
erosion_kge_bestObj = [erosion_nse(erosion_kge_sortindx(1,2),:) erosion_pbias(erosion_kge_sortindx(1,2),:) erosion_kge(erosion_kge_sortindx(1,2),:)];
% MULTI
% Normalise the objectives
erosion_NSEnorm   = (erosion_NSEmod-min(erosion_NSEmod))/(max(erosion_NSEmod)-min(erosion_NSEmod));
erosion_PBIASnorm = (erosion_PBIASmod-min(erosion_PBIASmod))/(max(erosion_PBIASmod)-min(erosion_PBIASmod));
erosion_KGEnorm = (erosion_KGEmod-min(erosion_KGEmod))/(max(erosion_KGEmod)-min(erosion_KGEmod));
erosion_multi = sqrt(erosion_NSEnorm.^2 + erosion_PBIASnorm.^2 + erosion_KGEnorm.^2);
% Best fit using Multi criteria
erosion_multi_sortindx = [];
[erosion_multi_sortindx(:,1),erosion_multi_sortindx(:,2)] = sort(erosion_multi);
erosion_multi_sortset = XY(erosion_multi_sortindx(:,2),:);
erosion_multi_bestObj = [erosion_nse(erosion_multi_sortindx(1,2),:) erosion_pbias(erosion_multi_sortindx(1,2),:) erosion_kge(erosion_multi_sortindx(1,2),:)];

% Simulations that comply with Moriasi et al. (2007) recommendations
erosion_M2007 = and(erosion_NSEmod<0.5,erosion_PBIASmod<55); % location
erosion_M2007n = 100*sum(erosion_M2007)/N3; % percentage

disp(' ')
%% 13. RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS
disp('RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS')

% RSA (find behavioural parameterizations):
erosion_threshold = [1-erosion_nse_thres erosion_pbias_thres 1-erosion_kge_thres];
[erosion_mvd_rsa,erosion_idxb_rsa] = RSA_indices_thres(XY,[erosion_NSEmod erosion_PBIASmod erosion_KGEmod],erosion_threshold);
erosion_Ybeh = XY(erosion_idxb_rsa,:);

erosion_top5n = 100*sum(erosion_idxb_rsa)/N3; % percentage of top 5% behavioural

% Highlight the behavioural parameterizations in the scatter plots:
figure; scatter_plots(Y,erosion_nse,[],'NSE',Y_Labels,erosion_idxb_rsa);
for i=M+1:M+E
    subplot(1,E,i-M)
    plot(erosion_nse_sortset(1,i),erosion_nse_bestObj(1),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_pbias_sortset(1,i),erosion_pbias_bestObj(1),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_kge_sortset(1,i),erosion_kge_bestObj(1),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_multi_sortset(1,i),erosion_multi_bestObj(1),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('Location','Southeast')

figure; scatter_plots(Y,erosion_pbias,[],'PBIAS',Y_Labels,erosion_idxb_rsa);
for i=M+1:M+E
    subplot(1,E,i-M)
    plot(erosion_nse_sortset(1,i),erosion_nse_bestObj(2),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_pbias_sortset(1,i),erosion_pbias_bestObj(2),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_kge_sortset(1,i),erosion_kge_bestObj(2),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_multi_sortset(1,i),erosion_multi_bestObj(2),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('Location','Southeast')

figure; scatter_plots(Y,erosion_kge,[],'KGE',Y_Labels,erosion_idxb_rsa);
for i=M+1:M+E
    subplot(1,E,i-M)
    plot(erosion_nse_sortset(1,i),erosion_nse_bestObj(3),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_pbias_sortset(1,i),erosion_pbias_bestObj(3),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_kge_sortset(1,i),erosion_kge_bestObj(3),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
    plot(erosion_multi_sortset(1,i),erosion_multi_bestObj(3),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1])
end
legend('Parameter population','Better than threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('Location','Southeast')

% Plot parameter CDFs:
RSA_plot_thres(Y,erosion_idxb_rsa,[],Y_Labels,{'behav','non-behav'}); % add legend

disp(' ')
%% 14. SENSITIVITY INDICES (MVD AND SPREAD) USING NSE, PBIAS, AND KGE
disp('SENSITIVITY INDICES (MVD AND SPREAD) USING NSE, PBIAS, AND KGE')

% mvd: maximum vertical distance between parameters CDFs

% Assess robustness by bootstrapping:
Nboot = 100;
[erosion_mvd,erosion_idxb,erosion_mvd_lb,erosion_mvd_ub] = RSA_indices_thres(Y,[erosion_NSEmod erosion_PBIASmod erosion_KGEmod],erosion_threshold,1,Nboot);
[erosion_spread,~,erosion_spread_lb,erosion_spread_ub] = RSA_indices_thres(Y,[erosion_NSEmod erosion_PBIASmod erosion_KGEmod],erosion_threshold,2,Nboot);
% Plot results:
figure; boxplot1(erosion_mvd,Y_Labels,'mvd',erosion_mvd_lb,erosion_mvd_ub)
figure; boxplot1(erosion_spread,Y_Labels,'spread',erosion_spread_lb,erosion_spread_ub)
set(gca,'YLimMode','auto')

disp(' ')
%% 15. PLOT ALL (HYDRO+EROSION) SENSITIVITY INDICES
disp('PLOT ALL (HYDRO+EROSION) SENSITIVITY INDICES')

% mvd: maximum vertical distance between parameters CDFs
% Plot results:
figure; boxplot1([hydro_multi_mvd,erosion_mvd],XY_Labels,'MVD',[hydro_multi_mvd_lb,erosion_mvd_lb],[hydro_multi_mvd_ub,erosion_mvd_ub])
grid on

disp(' ')
%% 16. COMPARING MODEL IMPLEMENTATIONS USING MANUAL OPTIMA
disp('COMPARING MODEL IMPLEMENTATIONS USING MANUAL OPTIMA')

% Manual optimum parameter sets
% Hydro:            A1  A2  A3 Sm1 Sm2 Sm3  BSm t(1/2) tau*
  manual_X(1,:) = [  0  20  60   0 150 250   70    70   20]; % Engda et al., 2012; Legesse et al. 2009
  manual_X(2,:) = [  2  14  50 200  10 100  100    70   10]; % Tilahun et al., 2013a; Tilahun et al., 2013b
  manual_X(3,:) = [  2  10  47 200  10  65  100    75   10]; % Guzman et al., 2017
  manual_X(4,:) = hydro_multi_sortset(1,:); % Best random set in hydro
% manual_X(5,:) = [ 25  40  15 400 100 700 4000 1/0.1 10*T]; % Collick et al., 2009

% Erosion:         at1 at2 as1 as2
  manual_Y(1,:) = [0.2 3.4 0.2 3.4]; % Tilahun et al. (2013a)
  manual_Y(2,:) = [  4   4   3   3]; % Tilahun et al. (2013b)
  manual_Y(3,:) = [  6   5 5.5 4.5]; % Guzman et al. (2017)
  manual_Y(4,:) = erosion_multi_sortset(1,M+1:M+E); % Best random set in erosion

% Run the model and compute selected model output:
manual_Qsim = model_evaluation(myfun,manual_X,rain,evap);
% Run the erosion model:
myfun3 = 'ethiopian_erosion';
[manual_Csim(1,:)] = model_evaluation(myfun3,manual_Y(1,:),0,[],manual_X(2,:),rain,evap); % Tilahun et al., 2013a
[manual_Csim(2,:)] = model_evaluation(myfun3,manual_Y(2,:),H2,[],manual_X(2,:),rain,evap); % Tilahun et al., 2013b
[manual_Csim(3,:)] = model_evaluation(myfun3,manual_Y(3,:),H,[],manual_X(3,:),rain,evap); % Guzman et al., 2017
[manual_Csim(4,:)] = model_evaluation(myfun2,erosion_multi_sortset(1,:),H,rain,evap); % Best random set in erosion
% [manual_Csim(5,:)] = model_evaluation(myfun2,manual_Y(4,:),H,[],manual_X(3,:),rain,evap); % Best random set + Guzman et al., 2017

% Non-positive sediment values are discarded
% sedm(~(sedm>0)) = nan;
% Only consider those values when sediment observations are present
manual_Csim(:,~(sedm>0)) = 0;

% Compute performance metrics
warmup = 90 ; % warmup period to be discarded
% Hydro:
manual_hydro_nse = NSE(manual_Qsim(:,warmup+1:end),flow(warmup+1:end)');
manual_hydro_pbias = PBIAS(manual_Qsim(:,warmup+1:end),flow(warmup+1:end)');
manual_hydro_kge = KGE(manual_Qsim(:,warmup+1:end),flow(warmup+1:end)');
% Erosion:
manual_erosion_nse = NSE(manual_Csim(:,warmup+1:end),sedm(warmup+1:end)');
manual_erosion_pbias = PBIAS(manual_Csim(:,warmup+1:end),sedm(warmup+1:end)');
manual_erosion_kge = KGE(manual_Csim(:,warmup+1:end),sedm(warmup+1:end)');

% % Plot one year per set from 1990 to 1993
% manual_Qsim(1,366:end) = NaN;
% manual_Qsim(2,[1:365,2*365+1:end]) = NaN;
% manual_Qsim(3,[1:2*365+1,T-365+1:end]) = NaN;
% manual_Qsim(4,1:T-365) = NaN;
% manual_Csim(1,366:end) = NaN;
% manual_Csim(2,[1:365,2*365+1:end]) = NaN;
% manual_Csim(3,[1:2*365+1,T-365+1:end]) = NaN;
% manual_Csim(4,1:T-365) = NaN;
% H(1:2*365+1) = NaN;
% H2([1:365,2*365+1:end]) = NaN;

% Plot simulated time series:
figure
% subplot(611)
subplot(6,1,1)
bar(rain,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
% Customise the plot
title('Rainfall')
ylabel('[mm d^{-1}]')
set(gca,'XLim',[0 T+1],'XTickLabel',[])
grid on

% subplot(612)
subplot(6,1,2)
plot(evap,'Color',[0.8 0.8 0.8],'LineWidth',2)
% Customise the plot
title('Potential evapotranspiration')
ylabel('[mm d^{-1}]')
set(gca,'XLim',[0 T+1],'XTickLabel',[])
grid on

% subplot(6,1,3:6)
subplot(6,1,3:4)
% add flow observations as black circles
a = plot(flow,'LineWidth',2,'Color',[0.8 0.8 0.8]);
hold on
% add manual optima
b = plot(manual_Qsim(1,:)','LineWidth',1,'LineStyle','--','Color',0.5*[1 1 1]);
c = plot(manual_Qsim(2,:)','LineWidth',1,'LineStyle','-.','Color',0.313725501298904*[1 1 1]);
d = plot(manual_Qsim(3,:)','LineWidth',1,'LineStyle','-','Color',[0 0 0]);
e = plot(manual_Qsim(4,:)','LineWidth',2,'LineStyle',':','Color',[0 0 0]);
% Customise the plot
title('Flow')
set(gca,'XLim',[0 T+1],'YLim',[-2.5 30],'XTickLabel',[])
ylabel('[mm d^{-1}]'); % xlabel('time [days]'); 
box on
grid on
legend('location','north')
legend([a b c d e],'Observed','Engda et al., 2012','Tilahun et al., 2013a, 2013b','Guzman et al., 2017b','Best sampled set')
legend('boxoff')
clear a b c d e

% Plot the erosion outputs
subplot(6,1,5:6)
% a = plot(flow,'Color',[0 0.447058826684952 0.74117648601532]);
a = plot(sedm,'Color',[0.8 0.8 0.8],'LineWidth',2);
hold on
% b = plot(erosion_Q,'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',1.5);
b = plot(manual_Csim(1,:)','LineWidth',1,'LineStyle','--','Color',0.313725501298904*[1 1 1]);
c = plot(manual_Csim(2,:)','LineWidth',1,'LineStyle','-.','Color',0.313725501298904*[1 1 1]);
d = plot(manual_Csim(3,:)','LineWidth',1,'LineStyle','-','Color',[0 0 0]);
e = plot(manual_Csim(4,:)','LineWidth',2,'LineStyle',':','Color',[0 0 0]);
% Fraction H
f = area(-5*H,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
plot(-5*H2,'-.k');
% Customise the plot
title('Sediment concentration')
set(gca,'XLim',[0 T+1],'YLim',[-6 40])
set(gca,'XTick',(0:50:T+1),'XTickLabels',datestr((0:50:T+1)+datenum('01-Jan-1990'),'dd/mm/yyyy'))
xlabel('time [day]'); ylabel('Concentration [g L^{-1}]')
box on
grid on
legend('location','north')
legend([a b c d e f],'Observed','Tilahun et al., 2013a','Tilahun et al., 2013b','Guzman et al., 2017b','Best sampled set','Fraction H (inverted)')
legend('boxoff')
clear a b c d e f g

disp(' ')
%% 17. PLOT BEHAVIOURAL PARAMETERISATIONS
disp('PLOT BEHAVIOURAL PARAMETERISATIONS')

figure(133)
% NSE CDF
subplot(4,3,1)
hydro_NSEyi = sort(hydro_NSEmod);
hydro_NSEFi = empiricalcdf(hydro_NSEmod,hydro_NSEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(hydro_NSEyi,hydro_NSEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_NSEyi(hydro_NSEyi<(1-hydro_nse_thres)),hydro_NSEFi(hydro_NSEyi<(1-hydro_nse_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[hydro_nse_thres hydro_nse_thres],[0 1],'--k');
plot(1-hydro_nse_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_nse_bestObj(1))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_pbias_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_pbias_bestObj(1))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_kge_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_kge_bestObj(1))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_multi_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_multi_bestObj(1))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-NSE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')

% PBIAS CDF
subplot(4,3,2)
hydro_PBIASyi = sort(hydro_PBIASmod);
hydro_PBIASFi = empiricalcdf(hydro_PBIASmod,hydro_PBIASyi);
% plot_cdf(PBIASmod,'abs(PBIAS) [%]')
% hold on
a = plot(hydro_PBIASyi,hydro_PBIASFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_PBIASyi(hydro_PBIASyi<hydro_pbias_thres),hydro_PBIASFi(hydro_PBIASyi<hydro_pbias_thres),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot([hydro_pbias_thres hydro_pbias_thres],[0 1],'--k');
plot(abs(hydro_nse_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_nse_bestObj(2))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_pbias_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_pbias_bestObj(2))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_kge_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_kge_bestObj(2))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_multi_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_multi_bestObj(2))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(0:10:50))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('abs(PBIAS) [%]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')
title('(a) Hydrological module')

% KGE CDF
subplot(4,3,3)
hydro_KGEyi = sort(hydro_KGEmod);
hydro_KGEFi = empiricalcdf(hydro_KGEmod,hydro_KGEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(hydro_KGEyi,hydro_KGEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_KGEyi(hydro_KGEyi<(1-hydro_kge_thres)),hydro_KGEFi(hydro_KGEyi<(1-hydro_kge_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[hydro_kge_thres hydro_kge_thres],[0 1],'--k');
d = plot(1-hydro_nse_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_nse_bestObj(3))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(1-hydro_pbias_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_pbias_bestObj(3))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(1-hydro_kge_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_kge_bestObj(3))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(1-hydro_multi_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_multi_bestObj(3))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-KGE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southwest')
legend([a b c d e f g],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('boxoff')
clear a b c d e f g
legend('off')

% NSE CDF
subplot(4,3,7)
erosion_NSEyi = sort(erosion_NSEmod);
erosion_NSEFi = empiricalcdf(erosion_NSEmod,erosion_NSEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(erosion_NSEyi,erosion_NSEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_NSEyi(erosion_NSEyi<(1-erosion_nse_thres)),erosion_NSEFi(erosion_NSEyi<(1-erosion_nse_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[erosion_nse_thres erosion_nse_thres],[0 1],'--k');
plot(1-erosion_nse_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_nse_bestObj(1))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_pbias_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_pbias_bestObj(1))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_kge_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_kge_bestObj(1))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_multi_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_multi_bestObj(1))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-NSE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')

% PBIAS CDF
subplot(4,3,8)
erosion_PBIASyi = sort(erosion_PBIASmod);
erosion_PBIASFi = empiricalcdf(erosion_PBIASmod,erosion_PBIASyi);
% plot_cdf(PBIASmod,'abs(PBIAS) [%]')
% hold on
a = plot(erosion_PBIASyi,erosion_PBIASFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_PBIASyi(erosion_PBIASyi<erosion_pbias_thres),erosion_PBIASFi(erosion_PBIASyi<erosion_pbias_thres),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot([erosion_pbias_thres erosion_pbias_thres],[0 1],'--k');
plot(abs(erosion_nse_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_nse_bestObj(2))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_pbias_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_pbias_bestObj(2))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_kge_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_kge_bestObj(2))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_multi_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_multi_bestObj(2))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(0:10:50))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('abs(PBIAS) [%]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')
title('(b) Erosion module')

% KGE CDF
subplot(4,3,9)
erosion_KGEyi = sort(erosion_KGEmod);
erosion_KGEFi = empiricalcdf(erosion_KGEmod,erosion_KGEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(erosion_KGEyi,erosion_KGEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_KGEyi(erosion_KGEyi<(1-erosion_kge_thres)),erosion_KGEFi(erosion_KGEyi<(1-erosion_kge_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[erosion_kge_thres erosion_kge_thres],[0 1],'--k');
d = plot(1-erosion_nse_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_nse_bestObj(3))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(1-erosion_pbias_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_pbias_bestObj(3))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(1-erosion_kge_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_kge_bestObj(3))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(1-erosion_multi_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_multi_bestObj(3))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-KGE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southwest')
legend([a b c d e f g],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('boxoff')
clear a b c d e f g

disp(' ')
%% 18. COORDINATE PLOTS OF BEHAVIOURAL SETS
disp('COORDINATE PLOTS OF BEHAVIOURAL SETS')

figure(133)
% Parallel coordinate plot:
subplot(4,3,4:5)
parcoor(X,X_Labels,[],hydro_multi_idxb_rsa,12)
hold on
multi_bestnorm_nse = (hydro_nse_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_pbias = (hydro_pbias_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_kge = (hydro_kge_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_multi = (hydro_multi_sortset(1,:)-min(X))./(max(X)-min(X));
plot(multi_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

% Parallel coordinate plot:
subplot(4,3,10:12)
parcoor(XY,XY_Labels,[],erosion_idxb_rsa,12)
hold on
erosion_bestnorm_nse = (erosion_nse_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_pbias = (erosion_pbias_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_kge = (erosion_kge_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_multi = (erosion_multi_sortset(1,:)-min(XY))./(max(XY)-min(XY));
plot(erosion_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

disp(' ')
%% 19. SCATTER PLOT BETWEEN NSE AND PBIAS
disp('SCATTER PLOT BETWEEN NSE AND PBIAS')

figure
subplot(2,1,1)
% Plot all simulations:
a = plot(hydro_nse,hydro_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
b = plot(hydro_nse(hydro_multi_idxb_rsa),hydro_pbias(hydro_multi_idxb_rsa),'ok','MarkerSize',4,'MarkerFaceColor','k');
% Threshold
c = plot([hydro_nse_thres hydro_nse_thres],hydro_pbias_thres*[-1 1],'--k');
plot([hydro_nse_thres 1],[hydro_pbias_thres hydro_pbias_thres],'--k');
plot([hydro_nse_thres 1],-1*[hydro_pbias_thres hydro_pbias_thres],'--k');
% Plot best predictions:
d = plot(hydro_nse_bestObj(1),hydro_nse_bestObj(2),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(hydro_pbias_bestObj(1),hydro_pbias_bestObj(2),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(hydro_kge_bestObj(1),hydro_kge_bestObj(2),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(hydro_multi_bestObj(1),hydro_multi_bestObj(2),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_hydro_nse(1),manual_hydro_pbias(1),'xk','MarkerSize',10);
text(manual_hydro_nse(1),manual_hydro_pbias(1),'  Engda et al. (2012)')
plot(manual_hydro_nse(2),manual_hydro_pbias(2),'xk','MarkerSize',10)
text(manual_hydro_nse(2),manual_hydro_pbias(2),'  Tilahun et al. (2013a, 2013b)')
plot(manual_hydro_nse(3),manual_hydro_pbias(3),'xk','MarkerSize',10)
text(manual_hydro_nse(3),manual_hydro_pbias(3),'  Guzman et al. (2017b)')
% Customise plot
set(gca,'XLim',[0.5 0.9],'YLim',[-40 40],'XTick',(0:0.05:1),'YTick',(-50:10:50))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(a) Hydrological module')
box on
grid on
grid minor
legend('location','northwest')
legend([a b c d e f g h],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima')
legend('boxoff')
clear a b c d e f g h
legend('off')

subplot(2,1,2)
% Plot all simulations:
a = plot(erosion_nse,erosion_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
b = plot(erosion_nse(erosion_idxb_rsa),erosion_pbias(erosion_idxb_rsa),'ok','MarkerSize',4,'MarkerFaceColor','k');
% Threshold
c = plot([erosion_nse_thres erosion_nse_thres],erosion_pbias_thres*[-1 1],'--k');
plot([erosion_nse_thres 1],[erosion_pbias_thres erosion_pbias_thres],'--k');
plot([erosion_nse_thres 1],-1*[erosion_pbias_thres erosion_pbias_thres],'--k');
% Plot best prediction:
d = plot(erosion_nse(erosion_nse_sortindx(1,2)),erosion_pbias(erosion_nse_sortindx(1,2)),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(erosion_nse(erosion_pbias_sortindx(1,2)),erosion_pbias(erosion_pbias_sortindx(1,2)),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(erosion_nse(erosion_kge_sortindx(1,2)),erosion_pbias(erosion_kge_sortindx(1,2)),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(erosion_nse(erosion_multi_sortindx(1,2)),erosion_pbias(erosion_multi_sortindx(1,2)),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_erosion_nse(1),manual_erosion_pbias(1),'xk','MarkerSize',10);
text(manual_erosion_nse(1),manual_erosion_pbias(1),'  Tilahun et al. (2013a)')
plot(manual_erosion_nse(2),manual_erosion_pbias(2),'xk','MarkerSize',10)
text(manual_erosion_nse(2),manual_erosion_pbias(2),'  Tilahun et al. (2013b)')
plot(manual_erosion_nse(3),manual_erosion_pbias(3),'xk','MarkerSize',10)
text(manual_erosion_nse(3),manual_erosion_pbias(3),'  Guzman et al. (2017b)')
% Customise plot
set(gca,'XLim',[0.4 0.8],'YLim',[-40 40],'XTick',(0:0.05:1),'YTick',(-50:10:50))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(b) Erosion module')
box on
grid on
grid minor
legend('location','northwest')
legend([a b c d e f g h],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima')
legend('boxoff')
clear a b c d e f g h

disp(' ')
%% 20. EXPERIMENT SETUP OF RANDOM COMBINATION OF HYDROLOGICAL AND EROSION MODULES
disp('RANDOM COMBINATION OF HYDROLOGICAL AND EROSION MODULES')

M  = 9 ; % number of uncertain parameters: 9
E  = 4; % number of erosion parameters: 4
ME  = M+E; % number of uncertain parameters: 13
DistrFun  = 'unif'  ; % Parameter distribution
% Parameter ranges:
xmin = [  0   5  20  50   5  50  50    20    1] ; % minimum values
xmax = [ 20  40  80 250 150 300 200   100   50] ; % maximum values
ymin = [       0         0      0     0] ; % minimum values
ymax = [       5         5     10    10] ; % maximum values
xymin = [xmin ymin] ; % minimum values
xymax = [xmax ymax] ; % maximum values

% Parameter names
X_Labels = {'A_1','A_2','A_3','S_{max_1}','S_{max_2}','S_{max_3}','BS_{max}','t_{1/2}','\tau^{*}'} ;
Y_Labels = {'a_{t_1}','a_{t_2}','a_{s_1}','a_{s_2}'} ;
XY_Labels = [X_Labels Y_Labels] ;

disp(' ')
% SAMPLING INPUT SPACE
disp('SAMPLING INPUT SPACE')

% Sample the input space using the 'AAT_sampling' function
% (see help AAT_sampling.m for more details)
SampStrategy = 'lhs' ; % sampling strategy (other options is 'rsu')
% Latin Hypercube Sampling
N = 10000 ; % sample size
DistrPar = cell(ME,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:ME; DistrPar{i} = [xymin(i) xymax(i)] ; end % this is to resort to
% the format required by AAT_sampling function.
% Perform sampling:
XY2 = AAT_sampling(SampStrategy,ME,DistrFun,DistrPar,N);
% Normalise the areas
% counter = 0; % Just for information
for i = 1:N
    if sum(XY2(i,1:3)) > 100
        % counter = counter + 1;
        XY2(i,1:3) = 100*XY2(i,1:3)/sum(XY2(i,1:3));
    end
end
% Adjust the erosion parameter values: at = (at-as)+as
XY2(:,M+1) = sum(XY2(:,[M+1 M+3]),2);
XY2(:,M+2) = sum(XY2(:,[M+2 M+4]),2);

% Plot results of the random sampling
figure; parcoor(XY2,XY_Labels); % Parallel coordinate plots

disp(' ')
%% 21. MODEL EVALUATION
disp('MODEL EVALUATION OF RANDOM COMBINATION')

% Erosion model:
myfun2 = 'ethiopian_ped';

% Run the model and compute model output at sampled parameter sets:
[C2sim,QC2sim,C2FLUXES] = model_evaluation(myfun2,XY2,H,rain,evap);

% Only consider those values when sediment observations are present
C2sim(:,~(sedm>0)) = 0;

% Extract fluxes and states:
Q2sim = nan(size(C2sim));
for i = 1:N
    Q2sim(i,:)     = QC2sim{i}(:,1)';
end

disp(' ')
%% 22. MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS
disp('MULTIOBJECTIVE SETUP AND REGIONAL SENSITIVITY ANALYSIS')

% Statistics of the erosion module
rand_erosion_nse = NSE(C2sim(:,warmup+1:end),sedm(warmup+1:end)');
rand_erosion_pbias = PBIAS(C2sim(:,warmup+1:end),sedm(warmup+1:end)');
rand_erosion_kge = KGE(C2sim(:,warmup+1:end),sedm(warmup+1:end)');

% Transformed objectives for minimisation
rand_erosion_NSEmod = 1 - (rand_erosion_nse); % 1-NSE to minimise
rand_erosion_PBIASmod = abs(rand_erosion_pbias); % abs(PBIAS) to minimise
rand_erosion_KGEmod = 1 - (rand_erosion_kge); % 1-KGE to minimise

% Set user-defined output thresholds
% rand_erosion_nse_thres = 0.7; % behavioural threshold for NSE
% rand_erosion_pbias_thres = 5; % behavioural threshold for PBIAS
% rand_erosion_kge_thres = 0.8; % % behavioural threshold for KGE, or same as NSE

% OR: Define the thresholds from the best values of each Obj
rand_erosion_nse_sortindx = []; rand_erosion_pbias_sortindx = []; rand_erosion_kge_sortindx = [];
[rand_erosion_nse_sortindx(:,1),rand_erosion_nse_sortindx(:,2)] = sort(rand_erosion_NSEmod);
[rand_erosion_pbias_sortindx(:,1),rand_erosion_pbias_sortindx(:,2)] = sort(rand_erosion_PBIASmod);
[rand_erosion_kge_sortindx(:,1),rand_erosion_kge_sortindx(:,2)] = sort(rand_erosion_KGEmod);
% Extract thresholds
rand_erosion_nse_thres = 1-rand_erosion_nse_sortindx(ceil(N*0.05));
rand_erosion_pbias_thres = abs(rand_erosion_pbias_sortindx(ceil(N*0.05)));
rand_erosion_kge_thres = 1-rand_erosion_kge_sortindx(ceil(N*0.05));

% Best fit using criteria:
% NSE
rand_erosion_nse_sortset = XY2(rand_erosion_nse_sortindx(:,2),:);
rand_erosion_nse_bestObj = [rand_erosion_nse(rand_erosion_nse_sortindx(1,2),:) rand_erosion_pbias(rand_erosion_nse_sortindx(1,2),:) rand_erosion_kge(rand_erosion_nse_sortindx(1,2),:)];
% PBIAS
rand_erosion_pbias_sortset = XY2(rand_erosion_pbias_sortindx(:,2),:);
rand_erosion_pbias_bestObj = [rand_erosion_nse(rand_erosion_pbias_sortindx(1,2),:) rand_erosion_pbias(rand_erosion_pbias_sortindx(1,2),:) rand_erosion_kge(rand_erosion_pbias_sortindx(1,2),:)];
% KGE
rand_erosion_kge_sortset = XY2(rand_erosion_kge_sortindx(:,2),:);
rand_erosion_kge_bestObj = [rand_erosion_nse(rand_erosion_kge_sortindx(1,2),:) rand_erosion_pbias(rand_erosion_kge_sortindx(1,2),:) rand_erosion_kge(rand_erosion_kge_sortindx(1,2),:)];
% MULTI
% Normalise the objectives
rand_erosion_NSEnorm   = (rand_erosion_NSEmod-min(rand_erosion_NSEmod))/(max(rand_erosion_NSEmod)-min(rand_erosion_NSEmod));
rand_erosion_PBIASnorm = (rand_erosion_PBIASmod-min(rand_erosion_PBIASmod))/(max(rand_erosion_PBIASmod)-min(rand_erosion_PBIASmod));
rand_erosion_KGEnorm = (rand_erosion_KGEmod-min(rand_erosion_KGEmod))/(max(rand_erosion_KGEmod)-min(rand_erosion_KGEmod));
rand_erosion_multi = sqrt(rand_erosion_NSEnorm.^2 + rand_erosion_PBIASnorm.^2 + rand_erosion_KGEnorm.^2);
% Best fit using Multi criteria
rand_erosion_multi_sortindx = [];
[rand_erosion_multi_sortindx(:,1),rand_erosion_multi_sortindx(:,2)] = sort(rand_erosion_multi);
rand_erosion_multi_sortset = XY2(rand_erosion_multi_sortindx(:,2),:);
rand_erosion_multi_bestObj = [rand_erosion_nse(rand_erosion_multi_sortindx(1,2),:) rand_erosion_pbias(rand_erosion_multi_sortindx(1,2),:) rand_erosion_kge(rand_erosion_multi_sortindx(1,2),:)];

disp(' ')
%% 23. EVALUATION OF THE RANDOM HYDROLOGICAL MODULE
disp('EVALUATION OF THE RANDOM HYDROLOGICAL MODULE')

% Statistics of the hydrological module
rand_hydro_nse = NSE(Q2sim(rand_erosion_idxb_rsa,warmup+1:end),flow(warmup+1:end)');
rand_hydro_pbias = PBIAS(Q2sim(rand_erosion_idxb_rsa,warmup+1:end),flow(warmup+1:end)');
rand_hydro_kge = KGE(Q2sim(rand_erosion_idxb_rsa,warmup+1:end),flow(warmup+1:end)');
% Transformed objectives for minimisation
rand_hydro_NSEmod = 1 - (rand_hydro_nse); % 1-NSE to minimise
rand_hydro_PBIASmod = abs(rand_hydro_pbias); % abs(PBIAS) to minimise
rand_hydro_KGEmod = 1 - (rand_hydro_kge); % 1-KGE to minimise

% OR: Define the thresholds from the best values of each Obj
rand_hydro_nse_sortindx = []; rand_hydro_pbias_sortindx = []; rand_hydro_kge_sortindx = [];
[rand_hydro_nse_sortindx(:,1),rand_hydro_nse_sortindx(:,2)] = sort(rand_hydro_NSEmod);
[rand_hydro_pbias_sortindx(:,1),rand_hydro_pbias_sortindx(:,2)] = sort(rand_hydro_PBIASmod);
[rand_hydro_kge_sortindx(:,1),rand_hydro_kge_sortindx(:,2)] = sort(rand_hydro_KGEmod);

% Best fit using criteria:
% NSE
rand_hydro_nse_sortset = rand_erosion_Ybeh(rand_hydro_nse_sortindx(:,2),:);
rand_hydro_nse_bestObj = [rand_hydro_nse(rand_hydro_nse_sortindx(1,2),:) rand_hydro_pbias(rand_hydro_nse_sortindx(1,2),:) rand_hydro_kge(rand_hydro_nse_sortindx(1,2),:)];
% PBIAS
rand_hydro_pbias_sortset = rand_erosion_Ybeh(rand_hydro_pbias_sortindx(:,2),:);
rand_hydro_pbias_bestObj = [rand_hydro_nse(rand_hydro_pbias_sortindx(1,2),:) rand_hydro_pbias(rand_hydro_pbias_sortindx(1,2),:) rand_hydro_kge(rand_hydro_pbias_sortindx(1,2),:)];
% KGE
rand_hydro_kge_sortset = rand_erosion_Ybeh(rand_hydro_kge_sortindx(:,2),:);
rand_hydro_kge_bestObj = [rand_hydro_nse(rand_hydro_kge_sortindx(1,2),:) rand_hydro_pbias(rand_hydro_kge_sortindx(1,2),:) rand_hydro_kge(rand_hydro_kge_sortindx(1,2),:)];
% MULTI
% Normalise the objectives
rand_hydro_NSEnorm   = (rand_hydro_NSEmod-min(rand_hydro_NSEmod))/(max(rand_hydro_NSEmod)-min(rand_hydro_NSEmod));
rand_hydro_PBIASnorm = (rand_hydro_PBIASmod-min(rand_hydro_PBIASmod))/(max(rand_hydro_PBIASmod)-min(rand_hydro_PBIASmod));
rand_hydro_KGEnorm = (rand_hydro_KGEmod-min(rand_hydro_KGEmod))/(max(rand_hydro_KGEmod)-min(rand_hydro_KGEmod));
rand_hydro_multi = sqrt(rand_hydro_NSEnorm.^2 + rand_hydro_PBIASnorm.^2 + rand_hydro_KGEnorm.^2);
% Best fit using Multi criteria
rand_hydro_multi_sortindx = [];
[rand_hydro_multi_sortindx(:,1),rand_hydro_multi_sortindx(:,2)] = sort(rand_hydro_multi);
rand_hydro_multi_sortset = rand_erosion_Ybeh(rand_hydro_multi_sortindx(:,2),:);
rand_hydro_multi_bestObj = [rand_hydro_nse(rand_hydro_multi_sortindx(1,2),:) rand_hydro_pbias(rand_hydro_multi_sortindx(1,2),:) rand_hydro_kge(rand_hydro_multi_sortindx(1,2),:)];

% Statistics of the hydrological module
rand_hydro_full_nse = NSE(Q2sim(:,warmup+1:end),flow(warmup+1:end)');
rand_hydro_full_pbias = PBIAS(Q2sim(:,warmup+1:end),flow(warmup+1:end)');
rand_hydro_full_kge = KGE(Q2sim(:,warmup+1:end),flow(warmup+1:end)');
% Transformed objectives for minimisation
rand_hydro_full_NSEmod = 1 - (rand_hydro_full_nse); % 1-NSE to minimise
rand_hydro_full_PBIASmod = abs(rand_hydro_full_pbias); % abs(PBIAS) to minimise
rand_hydro_full_KGEmod = 1 - (rand_hydro_full_kge); % 1-KGE to minimise

% Set user-defined output thresholds
% rand_hydro_full_nse_thres = 0.7; % behavioural threshold for NSE
% rand_hydro_full_pbias_thres = 5; % behavioural threshold for PBIAS
% rand_hydro_full_kge_thres = 0.8; % % behavioural threshold for KGE, or same as NSE

% OR: Define the thresholds from the best values of each Obj
rand_hydro_full_nse_sortindx = []; rand_hydro_full_pbias_sortindx = []; rand_hydro_full_kge_sortindx = [];
[rand_hydro_full_nse_sortindx(:,1),rand_hydro_full_nse_sortindx(:,2)] = sort(rand_hydro_full_NSEmod);
[rand_hydro_full_pbias_sortindx(:,1),rand_hydro_full_pbias_sortindx(:,2)] = sort(rand_hydro_full_PBIASmod);
[rand_hydro_full_kge_sortindx(:,1),rand_hydro_full_kge_sortindx(:,2)] = sort(rand_hydro_full_KGEmod);
% Extract thresholds
rand_hydro_full_nse_thres = 1-rand_hydro_full_nse_sortindx(ceil(N*0.05));
rand_hydro_full_pbias_thres = abs(rand_hydro_full_pbias_sortindx(ceil(N*0.05)));
rand_hydro_full_kge_thres = 1-rand_hydro_full_kge_sortindx(ceil(N*0.05));

% Best fit using criteria:
% NSE
rand_hydro_full_nse_sortset = XY2(rand_hydro_full_nse_sortindx(:,2),:);
rand_hydro_full_nse_bestObj = [rand_hydro_full_nse(rand_hydro_full_nse_sortindx(1,2),:) rand_hydro_full_pbias(rand_hydro_full_nse_sortindx(1,2),:) rand_hydro_full_kge(rand_hydro_full_nse_sortindx(1,2),:)];
% PBIAS
rand_hydro_full_pbias_sortset = XY2(rand_hydro_full_pbias_sortindx(:,2),:);
rand_hydro_full_pbias_bestObj = [rand_hydro_full_nse(rand_hydro_full_pbias_sortindx(1,2),:) rand_hydro_full_pbias(rand_hydro_full_pbias_sortindx(1,2),:) rand_hydro_full_kge(rand_hydro_full_pbias_sortindx(1,2),:)];
% KGE
rand_hydro_full_kge_sortset = XY2(rand_hydro_full_kge_sortindx(:,2),:);
rand_hydro_full_kge_bestObj = [rand_hydro_full_nse(rand_hydro_full_kge_sortindx(1,2),:) rand_hydro_full_pbias(rand_hydro_full_kge_sortindx(1,2),:) rand_hydro_full_kge(rand_hydro_full_kge_sortindx(1,2),:)];
% MULTI
% Normalise the objectives
rand_hydro_full_NSEnorm   = (rand_hydro_full_NSEmod-min(rand_hydro_full_NSEmod))/(max(rand_hydro_full_NSEmod)-min(rand_hydro_full_NSEmod));
rand_hydro_full_PBIASnorm = (rand_hydro_full_PBIASmod-min(rand_hydro_full_PBIASmod))/(max(rand_hydro_full_PBIASmod)-min(rand_hydro_full_PBIASmod));
rand_hydro_full_KGEnorm = (rand_hydro_full_KGEmod-min(rand_hydro_full_KGEmod))/(max(rand_hydro_full_KGEmod)-min(rand_hydro_full_KGEmod));% rand_hydro_full_multi = sqrt(rand_hydro_full_NSEnorm.^2 + rand_hydro_full_PBIASnorm.^2 + rand_hydro_full_KGEnorm.^2);
rand_hydro_full_multi = sqrt(rand_hydro_full_NSEnorm.^2 + rand_hydro_full_PBIASnorm.^2 + rand_hydro_full_KGEnorm.^2);
% Best fit using Multi criteria
rand_hydro_full_multi_sortindx = [];
[rand_hydro_full_multi_sortindx(:,1),rand_hydro_full_multi_sortindx(:,2)] = sort(rand_hydro_full_multi);
rand_hydro_full_multi_sortset = XY2(rand_hydro_full_multi_sortindx(:,2),:);
rand_hydro_full_multi_bestObj = [rand_hydro_full_nse(rand_hydro_full_multi_sortindx(1,2),:) rand_hydro_full_pbias(rand_hydro_full_multi_sortindx(1,2),:) rand_hydro_full_kge(rand_hydro_full_multi_sortindx(1,2),:)];

% Simulations that comply with Moriasi et al. (2007) recommendations
rand_erosion_M2007 = and(rand_erosion_NSEmod<0.5,rand_erosion_PBIASmod<55); % location
rand_erosion_M2007n = 100*sum(rand_erosion_M2007)/N; % percentage
rand_hydro_M2007 = and(rand_hydro_NSEmod<0.5,rand_hydro_PBIASmod<25); % location
rand_hydro_M2007n = 100*sum(rand_hydro_M2007)/N; % percentage
rand_hydro_full_M2007 = and(rand_hydro_full_NSEmod<0.5,rand_hydro_full_PBIASmod<25); % location
rand_hydro_full_M2007n = 100*sum(rand_hydro_full_M2007)/N; % percentage

disp(' ')
%% 24. RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS
disp('RSA USING NSE, PBIAS, AND KGE AS OBJECTIVE FUNCTIONS')

% RSA (find behavioural parameterizations):
rand_erosion_threshold = [1-rand_erosion_nse_thres rand_erosion_pbias_thres 1-rand_erosion_kge_thres];
[rand_erosion_mvd_rsa,rand_erosion_idxb_rsa] = RSA_indices_thres(XY2,[rand_erosion_NSEmod rand_erosion_PBIASmod rand_erosion_KGEmod],rand_erosion_threshold);
rand_erosion_Ybeh = XY2(rand_erosion_idxb_rsa,:);

rand_erosion_top5n = 100*sum(rand_erosion_idxb_rsa)/N; % percentage of top 5% behavioural

% Find behavioural parameterizations in the hydrological module from the ones in the erosion module:
rand_hydro_top5 = and(and(rand_hydro_NSEmod<(1-rand_hydro_full_nse_thres),rand_hydro_PBIASmod<rand_hydro_full_pbias_thres),rand_hydro_KGEmod<(1-rand_hydro_full_kge_thres)); % location
rand_hydro_top5n = 100*sum(rand_hydro_top5)/N; % percentage of top 5% behavioural

disp(' ')
%% 25. SENSITIVITY INDICES (MVD AND SPREAD) USING NSE, PBIAS, AND KGE
disp('SENSITIVITY INDICES (MVD AND SPREAD) USING NSE AND PBIAS')

% RSA (find behavioural parameterizations):
rand_erosion_threshold = [1-rand_erosion_nse_thres rand_erosion_pbias_thres 1-rand_erosion_kge_thres];
% mvd: maximum vertical distance between parameters CDFs

% Assess robustness by bootstrapping:
Nboot = 100;
[rand_erosion_mvd,rand_erosion_idxb,rand_erosion_mvd_lb,rand_erosion_mvd_ub] = RSA_indices_thres(XY2,[rand_erosion_NSEmod rand_erosion_PBIASmod rand_erosion_KGEmod],rand_erosion_threshold,1,Nboot);
[rand_erosion_spread,~,rand_erosion_spread_lb,rand_erosion_spread_ub] = RSA_indices_thres(XY2,[rand_erosion_NSEmod rand_erosion_PBIASmod rand_erosion_KGEmod],rand_erosion_threshold,2,Nboot);

disp(' ')
%% 26. PLOT SENSITIVITY INDICES (MVD AND SPREAD) OF BEHAVIOURAL AND RANDOM SETUPS
disp('PLOT SENSITIVITY INDICES (MVD AND SPREAD) OF BEHAVIOURAL AND RANDOM SETUPS')

% Plot results:
figure
subplot(2,2,1)
boxplot1([hydro_multi_mvd,erosion_mvd],XY_Labels,'mvd',[hydro_multi_mvd_lb,erosion_mvd_lb],[hydro_multi_mvd_ub,erosion_mvd_ub],12)
grid on
title('(a) Sensitivity using behavioural hydrological sets')

subplot(2,2,3)
boxplot1([hydro_multi_spread,erosion_spread],XY_Labels,'spread',[hydro_multi_spread_lb,erosion_spread_lb],[hydro_multi_spread_ub,erosion_spread_ub],12)
grid on
set(gca,'YLimMode','auto')

subplot(2,2,2)
boxplot1(rand_erosion_mvd,XY_Labels,'mvd',rand_erosion_mvd_lb,rand_erosion_mvd_ub,12)
grid on
title('(b) Sensitivity using random hydrological sets')

subplot(2,2,4)
boxplot1(rand_erosion_spread,XY_Labels,'spread',rand_erosion_spread_lb,rand_erosion_spread_ub,12)
grid on
set(gca,'YLimMode','auto')

disp(' ')
%% 27. SCATTER PLOT BETWEEN NSE AND PBIAS FOR THE RANDOM MODEL
disp('SCATTER PLOT BETWEEN NSE AND PBIAS')

figure
% subplot(2,2,1)
subplot(7,7,[1:3 8:10 15:17])
% Plot all simulations:
a = plot(hydro_nse,hydro_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
% Threshold
c = plot([hydro_nse_thres hydro_nse_thres],hydro_pbias_thres*[-1 1],'--k');
plot([hydro_nse_thres 1],[hydro_pbias_thres hydro_pbias_thres],'--k');
plot([hydro_nse_thres 1],-1*[hydro_pbias_thres hydro_pbias_thres],'--k');
% Behavioural sets
b = plot(hydro_nse(hydro_multi_idxb_rsa),hydro_pbias(hydro_multi_idxb_rsa),'ok','MarkerSize',4,'MarkerFaceColor','k');
% Plotting behavioural erosion sets into the hydrolgical module
i = plot(2,0,'sk','MarkerFaceColor','k','MarkerSize',4);
% Plot best predictions:
d = plot(hydro_nse_bestObj(1),hydro_nse_bestObj(2),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(hydro_pbias_bestObj(1),hydro_pbias_bestObj(2),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(hydro_kge_bestObj(1),hydro_kge_bestObj(2),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(hydro_multi_bestObj(1),hydro_multi_bestObj(2),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_hydro_nse(1),manual_hydro_pbias(1),'xk','MarkerSize',10,'LineWidth',2);
text(manual_hydro_nse(1),manual_hydro_pbias(1),'  Engda et al. (2012)')
plot(manual_hydro_nse(2),manual_hydro_pbias(2),'xk','MarkerSize',10,'LineWidth',2)
text(manual_hydro_nse(2),manual_hydro_pbias(2),'  Tilahun et al. (2013a, 2013b)')
plot(manual_hydro_nse(3),manual_hydro_pbias(3),'xk','MarkerSize',10,'LineWidth',2)
text(manual_hydro_nse(3),manual_hydro_pbias(3),'  Guzman et al. (2017b)')
% Thresholds from Moriasi et al. (2007)
plot([0.5 0.5],25*[-1 1],'--k');
plot([0.5 1],[25 25],'--k');
plot([0.5 1],-1*[25 25],'--k');
% text(0.41,5,'  Moriasi et al. (2007) threshold')
% Customise plot
set(gca,'XLim',[0.4 0.9],'YLim',[-50 50],'XTick',(0:0.05:1),'YTick',(-100:20:100))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(a) Hydrological module')
box on
grid on
% grid minor
% Legend 
legend('location','northwest')
legend('Orientation','horizontal')
legx = legend([a c b i d e f g h],'Parameter population','Threshold','Behavioural hydrological sets','Behavioural erosion sets','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima');
% legend('boxoff')
lpox = get(legx,'Position');
lpox(2) = 0.49;
set(legx,'Position',lpox)
clear a b c d e f g h i lpox
% legend('off')

% subplot(2,2,3)
subplot(7,7,[29:31 36:38 43:45])
% Plot all simulations:
a = plot(rand_hydro_full_nse,rand_hydro_full_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
% Threshold
c = plot([rand_hydro_full_nse_thres rand_hydro_full_nse_thres],rand_hydro_full_pbias_thres*[-1 1],'--k');
plot([rand_hydro_full_nse_thres 1],[rand_hydro_full_pbias_thres rand_hydro_full_pbias_thres],'--k');
plot([rand_hydro_full_nse_thres 1],-1*[rand_hydro_full_pbias_thres rand_hydro_full_pbias_thres],'--k');
% Behavioural hydrological sets
% b = plot(rand_hydro_full_nse(rand_hydro_full_multi_idxb_rsa),rand_hydro_full_pbias(rand_hydro_full_multi_idxb_rsa),'ok','MarkerSize',4,'MarkerFaceColor','k');
% Plotting behavioural erosion into hydrolgical module
i = plot(rand_hydro_nse,rand_hydro_pbias,'sk','MarkerSize',4,'MarkerFaceColor','k');
% Plot best predictions:
d = plot(rand_hydro_full_nse_bestObj(1),rand_hydro_full_nse_bestObj(2),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(rand_hydro_full_pbias_bestObj(1),rand_hydro_full_pbias_bestObj(2),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(rand_hydro_full_kge_bestObj(1),rand_hydro_full_kge_bestObj(2),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(rand_hydro_full_multi_bestObj(1),rand_hydro_full_multi_bestObj(2),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_hydro_nse(1),manual_hydro_pbias(1),'xk','MarkerSize',10,'LineWidth',2);
text(manual_hydro_nse(1),manual_hydro_pbias(1),'  Engda et al. (2012)')
plot(manual_hydro_nse(2),manual_hydro_pbias(2),'xk','MarkerSize',10,'LineWidth',2)
text(manual_hydro_nse(2),manual_hydro_pbias(2),'  Tilahun et al. (2013a, 2013b)')
plot(manual_hydro_nse(3),manual_hydro_pbias(3),'xk','MarkerSize',10,'LineWidth',2)
text(manual_hydro_nse(3),manual_hydro_pbias(3),'  Guzman et al. (2017b)')
% Thresholds from Moriasi et al. (2007)
plot([0.5 0.5],25*[-1 1],'--k');
plot([0.5 1],[25 25],'--k');
plot([0.5 1],-1*[25 25],'--k');
% text(0.5,0,'  Moriasi et al. (2007) threshold')
% Customise plot
set(gca,'XLim',[0.4 0.9],'YLim',[-50 50],'XTick',(0:0.05:1),'YTick',(-100:20:100))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(d) Hydrological module using behavioural erosion sets')
box on
grid on
% grid minor
% Legend 
legend('location','northwest')
legend([a c i d e f g h],'Parameter population','Threshold','Behavioural erosion sets','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima')
legend('boxoff')
clear a b c d e f g h i
legend('off')

% subplot(2,2,2)
subplot(7,7,[5:7 12:14 19:21])
% Plot all simulations:
a = plot(erosion_nse,erosion_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
% Threshold
c = plot([erosion_nse_thres erosion_nse_thres],erosion_pbias_thres*[-1 1],'--k');
plot([erosion_nse_thres 1],[erosion_pbias_thres erosion_pbias_thres],'--k');
plot([erosion_nse_thres 1],-1*[erosion_pbias_thres erosion_pbias_thres],'--k');
% Behavioural erosion sets
b = plot(erosion_nse(erosion_idxb_rsa),erosion_pbias(erosion_idxb_rsa),'ok','MarkerSize',4,'MarkerFaceColor','k');
% Plot best prediction:
d = plot(erosion_nse(erosion_nse_sortindx(1,2)),erosion_pbias(erosion_nse_sortindx(1,2)),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(erosion_nse(erosion_pbias_sortindx(1,2)),erosion_pbias(erosion_pbias_sortindx(1,2)),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(erosion_nse(erosion_kge_sortindx(1,2)),erosion_pbias(erosion_kge_sortindx(1,2)),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(erosion_nse(erosion_multi_sortindx(1,2)),erosion_pbias(erosion_multi_sortindx(1,2)),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_erosion_nse(1),manual_erosion_pbias(1),'xk','MarkerSize',10,'LineWidth',2);
text(manual_erosion_nse(1),manual_erosion_pbias(1),'  Tilahun et al. (2013a)')
plot(manual_erosion_nse(2),manual_erosion_pbias(2),'xk','MarkerSize',10,'LineWidth',2)
text(manual_erosion_nse(2),manual_erosion_pbias(2),'  Tilahun et al. (2013b)')
plot(manual_erosion_nse(3),manual_erosion_pbias(3),'xk','MarkerSize',10,'LineWidth',2)
text(manual_erosion_nse(3),manual_erosion_pbias(3),'  Guzman et al. (2017b)')
% Thresholds from Moriasi et al. (2007)
plot([0.5 0.5],55*[-1 1],'--k');
plot([0.5 1],[55 55],'--k');
plot([0.5 1],-1*[55 55],'--k');
% text(0.5,0,'  Moriasi et al. (2007) threshold')
% Customise plot
set(gca,'XLim',[0.3 0.8],'YLim',[-50 50],'XTick',(0:0.05:1),'YTick',(-100:20:100))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(b) Erosion module using behavioural hydrological sets')
box on
grid on
% grid minor
legend('location','southeast')
legend([a c b d e f g h],'Parameter population','Threshold','Behavioural erosion sets','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima')
legend('boxoff')
clear a b c d e f g h
legend('off')

% subplot(2,2,4)
subplot(7,7,[33:35 40:42 47:49])
% Plot all simulations:
a = plot(rand_erosion_nse,rand_erosion_pbias,'o','MarkerSize',3,'Color',[.5 .5 .5]);
hold on
% Threshold
c = plot([rand_erosion_nse_thres rand_erosion_nse_thres],rand_erosion_pbias_thres*[-1 1],'--k');
plot([rand_erosion_nse_thres 1],[rand_erosion_pbias_thres rand_erosion_pbias_thres],'--k');
plot([rand_erosion_nse_thres 1],-1*[rand_erosion_pbias_thres rand_erosion_pbias_thres],'--k');
% text(0.5,0,'  Moriasi et al. (2007) threshold')
% Behavioural erosion parameters
i = plot(rand_erosion_nse(rand_erosion_idxb_rsa),rand_erosion_pbias(rand_erosion_idxb_rsa),'sk','MarkerSize',4,'MarkerFaceColor','k');
% Plot best prediction:
d = plot(rand_erosion_nse(rand_erosion_nse_sortindx(1,2)),rand_erosion_pbias(rand_erosion_nse_sortindx(1,2)),'vk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(rand_erosion_nse(rand_erosion_pbias_sortindx(1,2)),rand_erosion_pbias(rand_erosion_pbias_sortindx(1,2)),'^k','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(rand_erosion_nse(rand_erosion_kge_sortindx(1,2)),rand_erosion_pbias(rand_erosion_kge_sortindx(1,2)),'sk','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(rand_erosion_nse(rand_erosion_multi_sortindx(1,2)),rand_erosion_pbias(rand_erosion_multi_sortindx(1,2)),'ok','MarkerSize',10,'MarkerFaceColor',0.8*[1 1 1]);
% Manual optima
h = plot(manual_erosion_nse(1),manual_erosion_pbias(1),'xk','MarkerSize',10,'LineWidth',2);
text(manual_erosion_nse(1),manual_erosion_pbias(1),'  Tilahun et al. (2013a)')
plot(manual_erosion_nse(2),manual_erosion_pbias(2),'xk','MarkerSize',10,'LineWidth',2)
text(manual_erosion_nse(2),manual_erosion_pbias(2),'  Tilahun et al. (2013b)')
plot(manual_erosion_nse(3),manual_erosion_pbias(3),'xk','MarkerSize',10,'LineWidth',2)
text(manual_erosion_nse(3),manual_erosion_pbias(3),'  Guzman et al. (2017b)')
% Thresholds from Moriasi et al. (2007)
plot([0.5 0.5],55*[-1 1],'--k');
plot([0.5 1],[55 55],'--k');
plot([0.5 1],-1*[55 55],'--k');
% Customise plot
set(gca,'XLim',[0.3 0.8],'YLim',[-50 50],'XTick',(0:0.05:1),'YTick',(-100:20:100))
set(gca,'FontName','Helvetica','FontSize',12)
xlabel('NSE [-]'); ylabel('PBIAS [%]')
title('(c) Erosion module using random hydrological sets')
box on
grid on
% grid minor
legend('location','northeast')
legend([a c i d e f g h],'Parameter population','Threshold','Behavioural erosion sets','Best NSE set','Best PBIAS set','Best KGE set','Best set overall','Manual optima')
legend('boxoff')
clear a b c d e f g h i
legend('off')

subplot(7,7,11)
img_arrow_right = imread('arrow_right.jpg');
image(img_arrow_right);
apox = get(gca,'Position');
apox(1) = apox(1)+apox(3)/4; apox(3)= apox(3)/2;
set(gca,'Color','none','Xtick',[],'Ytick',[],'Position',apox)
clear apox img_arrow_right
axis off

subplot(7,7,39)
img_arrow_left = imread('arrow_left.jpg');
image(img_arrow_left);
apox = get(gca,'Position');
apox(1) = apox(1)+apox(3)/4; apox(3)= apox(3)/2;
set(gca,'Color','none','Xtick',[],'Ytick',[],'Position',apox)
clear apox img_arrow_left
axis off

disp(' ')
%% 28. PLOT BEHAVIOURAL PARAMETERISATIONS WITH THE RANDOM MODEL
disp('PLOT BEHAVIOURAL PARAMETERISATIONS')

figure(134)

% HYDROLOGICAL MODULE
% NSE CDF
subplot(6,3,1)
hydro_NSEyi = sort(hydro_NSEmod);
hydro_NSEFi = empiricalcdf(hydro_NSEmod,hydro_NSEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(hydro_NSEyi,hydro_NSEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_NSEyi(hydro_NSEyi<(1-hydro_nse_thres)),hydro_NSEFi(hydro_NSEyi<(1-hydro_nse_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[hydro_nse_thres hydro_nse_thres],[0 1],'--k');
plot(1-hydro_nse_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_nse_bestObj(1))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_pbias_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_pbias_bestObj(1))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_kge_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_kge_bestObj(1))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-hydro_multi_bestObj(1),hydro_NSEFi(hydro_NSEyi==(1-hydro_multi_bestObj(1))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-NSE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')

% PBIAS CDF
subplot(6,3,2)
hydro_PBIASyi = sort(hydro_PBIASmod);
hydro_PBIASFi = empiricalcdf(hydro_PBIASmod,hydro_PBIASyi);
% plot_cdf(PBIASmod,'abs(PBIAS) [%]')
% hold on
a = plot(hydro_PBIASyi,hydro_PBIASFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_PBIASyi(hydro_PBIASyi<hydro_pbias_thres),hydro_PBIASFi(hydro_PBIASyi<hydro_pbias_thres),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot([hydro_pbias_thres hydro_pbias_thres],[0 1],'--k');
plot(abs(hydro_nse_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_nse_bestObj(2))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_pbias_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_pbias_bestObj(2))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_kge_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_kge_bestObj(2))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(hydro_multi_bestObj(2)),hydro_PBIASFi(hydro_PBIASyi==abs(hydro_multi_bestObj(2))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(0:10:50))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('abs(PBIAS) [%]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')
title('(a) Hydrological module')

% KGE CDF
subplot(6,3,3)
hydro_KGEyi = sort(hydro_KGEmod);
hydro_KGEFi = empiricalcdf(hydro_KGEmod,hydro_KGEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(hydro_KGEyi,hydro_KGEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(hydro_KGEyi(hydro_KGEyi<(1-hydro_kge_thres)),hydro_KGEFi(hydro_KGEyi<(1-hydro_kge_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[hydro_kge_thres hydro_kge_thres],[0 1],'--k');
d = plot(1-hydro_nse_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_nse_bestObj(3))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(1-hydro_pbias_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_pbias_bestObj(3))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(1-hydro_kge_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_kge_bestObj(3))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(1-hydro_multi_bestObj(3),hydro_KGEFi(hydro_KGEyi==(1-hydro_multi_bestObj(3))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-KGE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southwest')
legend([a b c d e f g],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('boxoff')
clear a b c d e f g
legend('off')

% EROSION MODULE
% NSE CDF
subplot(6,3,7)
erosion_NSEyi = sort(erosion_NSEmod);
erosion_NSEFi = empiricalcdf(erosion_NSEmod,erosion_NSEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(erosion_NSEyi,erosion_NSEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_NSEyi(erosion_NSEyi<(1-erosion_nse_thres)),erosion_NSEFi(erosion_NSEyi<(1-erosion_nse_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[erosion_nse_thres erosion_nse_thres],[0 1],'--k');
plot(1-erosion_nse_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_nse_bestObj(1))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_pbias_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_pbias_bestObj(1))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_kge_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_kge_bestObj(1))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-erosion_multi_bestObj(1),erosion_NSEFi(erosion_NSEyi==(1-erosion_multi_bestObj(1))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
erosion_NSEYLim = get(gca,'XLim');
set(gca,'XLim',erosion_NSEYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-NSE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')

% PBIAS CDF
subplot(6,3,8)
erosion_PBIASyi = sort(erosion_PBIASmod);
erosion_PBIASFi = empiricalcdf(erosion_PBIASmod,erosion_PBIASyi);
% plot_cdf(PBIASmod,'abs(PBIAS) [%]')
% hold on
a = plot(erosion_PBIASyi,erosion_PBIASFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_PBIASyi(erosion_PBIASyi<erosion_pbias_thres),erosion_PBIASFi(erosion_PBIASyi<erosion_pbias_thres),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot([erosion_pbias_thres erosion_pbias_thres],[0 1],'--k');
plot(abs(erosion_nse_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_nse_bestObj(2))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_pbias_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_pbias_bestObj(2))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_kge_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_kge_bestObj(2))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(erosion_multi_bestObj(2)),erosion_PBIASFi(erosion_PBIASyi==abs(erosion_multi_bestObj(2))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(0:10:50))
erosion_PBIASYLim = get(gca,'XLim');
set(gca,'XLim',erosion_PBIASYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('abs(PBIAS) [%]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c
legend('off')
title('(b) Erosion module using behavioural hydrological sets')

% KGE CDF
subplot(6,3,9)
erosion_KGEyi = sort(erosion_KGEmod);
erosion_KGEFi = empiricalcdf(erosion_KGEmod,erosion_KGEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(erosion_KGEyi,erosion_KGEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(erosion_KGEyi(erosion_KGEyi<(1-erosion_kge_thres)),erosion_KGEFi(erosion_KGEyi<(1-erosion_kge_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[erosion_kge_thres erosion_kge_thres],[0 1],'--k');
d = plot(1-erosion_nse_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_nse_bestObj(3))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(1-erosion_pbias_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_pbias_bestObj(3))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(1-erosion_kge_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_kge_bestObj(3))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(1-erosion_multi_bestObj(3),erosion_KGEFi(erosion_KGEyi==(1-erosion_multi_bestObj(3))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
erosion_KGEYLim = get(gca,'XLim');
set(gca,'XLim',erosion_KGEYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-KGE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southwest')
legend([a b c d e f g],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('boxoff')
clear a b c d e f g
legend('off')

% RANDOM EROSION MODULE
% NSE CDF
subplot(6,3,13)
rand_erosion_NSEyi = sort(rand_erosion_NSEmod);
rand_erosion_NSEFi = empiricalcdf(rand_erosion_NSEmod,rand_erosion_NSEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(rand_erosion_NSEyi,rand_erosion_NSEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(rand_erosion_NSEyi(rand_erosion_NSEyi<(1-rand_erosion_nse_thres)),rand_erosion_NSEFi(rand_erosion_NSEyi<(1-rand_erosion_nse_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[rand_erosion_nse_thres rand_erosion_nse_thres],[0 1],'--k');
plot(1-rand_erosion_nse_bestObj(1),rand_erosion_NSEFi(rand_erosion_NSEyi==(1-rand_erosion_nse_bestObj(1))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-rand_erosion_pbias_bestObj(1),rand_erosion_NSEFi(rand_erosion_NSEyi==(1-rand_erosion_pbias_bestObj(1))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-rand_erosion_kge_bestObj(1),rand_erosion_NSEFi(rand_erosion_NSEyi==(1-rand_erosion_kge_bestObj(1))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(1-rand_erosion_multi_bestObj(1),rand_erosion_NSEFi(rand_erosion_NSEyi==(1-rand_erosion_multi_bestObj(1))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
set(gca,'XLim',erosion_NSEYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-NSE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c erosion_NSEYLim
legend('off')

% PBIAS CDF
subplot(6,3,14)
rand_erosion_PBIASyi = sort(rand_erosion_PBIASmod);
rand_erosion_PBIASFi = empiricalcdf(rand_erosion_PBIASmod,rand_erosion_PBIASyi);
% plot_cdf(PBIASmod,'abs(PBIAS) [%]')
% hold on
a = plot(rand_erosion_PBIASyi,rand_erosion_PBIASFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(rand_erosion_PBIASyi(rand_erosion_PBIASyi<rand_erosion_pbias_thres),rand_erosion_PBIASFi(rand_erosion_PBIASyi<rand_erosion_pbias_thres),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot([rand_erosion_pbias_thres rand_erosion_pbias_thres],[0 1],'--k');
plot(abs(rand_erosion_nse_bestObj(2)),rand_erosion_PBIASFi(rand_erosion_PBIASyi==abs(rand_erosion_nse_bestObj(2))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(rand_erosion_pbias_bestObj(2)),rand_erosion_PBIASFi(rand_erosion_PBIASyi==abs(rand_erosion_pbias_bestObj(2))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(rand_erosion_kge_bestObj(2)),rand_erosion_PBIASFi(rand_erosion_PBIASyi==abs(rand_erosion_kge_bestObj(2))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
plot(abs(rand_erosion_multi_bestObj(2)),rand_erosion_PBIASFi(rand_erosion_PBIASyi==abs(rand_erosion_multi_bestObj(2))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(0:10:50))
set(gca,'XLim',erosion_PBIASYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('abs(PBIAS) [%]'); ylabel('CDF')
grid on
grid minor
legend('location','southeast')
legend([a b c],'Parameter population','Better than threshold','Threshold')
clear a b c erosion_PBIASYLim
legend('off')
title('(c) Erosion module using random hydrological sets')

% KGE CDF
subplot(6,3,15)
rand_erosion_KGEyi = sort(rand_erosion_KGEmod);
rand_erosion_KGEFi = empiricalcdf(rand_erosion_KGEmod,rand_erosion_KGEyi);
% plot_cdf(NSEmod,'1-NSE [-]');
% hold on
a = plot(rand_erosion_KGEyi,rand_erosion_KGEFi,'o','MarkerSize',3,'Color',0.5*[1 1 1]);
hold on
b = plot(rand_erosion_KGEyi(rand_erosion_KGEyi<(1-rand_erosion_kge_thres)),rand_erosion_KGEFi(rand_erosion_KGEyi<(1-rand_erosion_kge_thres)),'ok','MarkerFaceColor','k','MarkerSize',4);
c = plot(1-[rand_erosion_kge_thres rand_erosion_kge_thres],[0 1],'--k');
d = plot(1-rand_erosion_nse_bestObj(3),rand_erosion_KGEFi(rand_erosion_KGEyi==(1-rand_erosion_nse_bestObj(3))),'vk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
e = plot(1-rand_erosion_pbias_bestObj(3),rand_erosion_KGEFi(rand_erosion_KGEyi==(1-rand_erosion_pbias_bestObj(3))),'^k','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
f = plot(1-rand_erosion_kge_bestObj(3),rand_erosion_KGEFi(rand_erosion_KGEyi==(1-rand_erosion_kge_bestObj(3))),'sk','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
g = plot(1-rand_erosion_multi_bestObj(3),rand_erosion_KGEFi(rand_erosion_KGEyi==(1-rand_erosion_multi_bestObj(3))),'ok','MarkerSize',6,'MarkerFaceColor',0.8*[1 1 1]);
% Customise the plot
set(gca,'YLim',[0 0.25],'YTick',(0:0.10:1),'XTick',(-1:0.10:1))
set(gca,'XLim',erosion_KGEYLim);
% set(gca,'FontName','Helvetica','FontSize',20)
xlabel('1-KGE [-]'); ylabel('CDF')
grid on
grid minor
legend('location','southwest')
legend([a b c d e f g],'Parameter population','Better than threshold','Threshold','Best NSE set','Best PBIAS set','Best KGE set','Best set overall')
legend('boxoff')
clear a b c d e f g erosion_KGEYLim

disp(' ')
%% 29. COORDINATE PLOTS OF BEHAVIOURAL SETS WITH THE RANDOM MODEL
disp('COORDINATE PLOTS OF BEHAVIOURAL SETS')

figure(134)
% Parallel coordinate plot:
subplot(6,3,4:5)
parcoor(X,X_Labels,[],hydro_multi_idxb_rsa,10)
hold on
multi_bestnorm_nse = (hydro_nse_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_pbias = (hydro_pbias_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_kge = (hydro_kge_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_multi = (hydro_multi_sortset(1,:)-min(X))./(max(X)-min(X));
plot(multi_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

% Parallel coordinate plot:
% subplot(6,3,10:12)
% parcoor(Y,Y_Labels,[],erosion_idxb_rsa,12)
% hold on
% erosion_bestnorm_nse = (erosion_nse_sortset(1,M+1:end)-min(Y))./(max(Y)-min(Y));
% erosion_bestnorm_pbias = (erosion_pbias_sortset(1,M+1:end)-min(Y))./(max(Y)-min(Y));
% erosion_bestnorm_kge = (erosion_kge_sortset(1,M+1:end)-min(Y))./(max(Y)-min(Y));
% erosion_bestnorm_multi = (erosion_multi_sortset(1,M+1:end)-min(Y))./(max(Y)-min(Y));
% plot(erosion_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
% plot(erosion_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
% plot(erosion_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
% plot(erosion_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
subplot(6,3,10:12)
parcoor(XY,XY_Labels,[],erosion_idxb_rsa,10)
hold on
erosion_bestnorm_nse = (erosion_nse_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_pbias = (erosion_pbias_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_kge = (erosion_kge_sortset(1,:)-min(XY))./(max(XY)-min(XY));
erosion_bestnorm_multi = (erosion_multi_sortset(1,:)-min(XY))./(max(XY)-min(XY));
plot(erosion_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(erosion_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

% Parallel coordinate plot:
subplot(6,3,16:18)
parcoor(XY2,XY_Labels,[],rand_erosion_idxb_rsa,10)
hold on
rand_erosion_bestnorm_nse = (rand_erosion_nse_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_pbias = (rand_erosion_pbias_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_kge = (rand_erosion_kge_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_multi = (rand_erosion_multi_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
plot(rand_erosion_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

disp(' ')
%% 30. SETUP PAWN FOR THE HYDROLOGICAL MODULE
disp('SETUP PAWN FOR THE HYDROLOGICAL MODULE')

% NU = 150 ; % number of samples to estimate unconditional CDF
NC = 1000 ; % number of samples to estimate conditional CDFs
n  = 10 ; % number of conditioning points

DistrPar = cell(M,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:M; DistrPar{i} = [xmin(i) xmax(i) ] ; end % this is to resort to
% the format required by AAT_sampling function.

% % Create input/output samples to estimate the unconditional output CDF:
% hydro_pawn_Xu = AAT_sampling('lhs',M,'unif',DistrPar,NU); % matrix (NU,M)
% % Normalise the areas
% % counter = 0; % Just for information
% for i = 1:NU
%     if sum(hydro_pawn_Xu(i,1:3)) > 100
%         % counter = counter + 1;
%         hydro_pawn_Xu(i,1:3) = 100*hydro_pawn_Xu(i,1:3)/sum(hydro_pawn_Xu(i,1:3));
%     end
% end
% hydro_pawn_Yu = model_evaluation(myfun,hydro_pawn_Xu,rain,evap)  ; % vector (1,M)

% Create input/output samples to estimate the conditional output CDFs:
[ hydro_pawn_XX, hydro_pawn_xc ] = pawn_sampling('lhs',M,'unif',DistrPar,n,NC);
% Normalise the areas
for j = 1:n
    for h = 1:M
        % counter = 0; % Just for information
        for i = 1:NC
            if sum(hydro_pawn_XX{h,j}(i,1:3)) > 100
                % counter = counter + 1;
                hydro_pawn_XX{h,j}(i,1:3) = 100*hydro_pawn_XX{h,j}(i,1:3)/sum(hydro_pawn_XX{h,j}(i,1:3));
            end
        end
    end
end
hydro_pawn_Qsim_YY = pawn_model_evaluation(myfun,hydro_pawn_XX,rain,evap) ;

disp(' ')
%% 31. EVALUATE PAWN STATISTICS FOR THE HYDROLOGICAL MODULE
disp('EVALUATE PAWN STATISTICS FOR THE HYDROLOGICAL MODULE')

% Statistics
hydro_pawn_nse = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_pbias = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_kge = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_NSEmod = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_PBIASmod = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_KGEmod = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_NSEnorm = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_PBIASnorm = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_KGEnorm = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_multi = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_maxQ = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_minQ = cell(size(hydro_pawn_Qsim_YY));
hydro_pawn_meanQ = cell(size(hydro_pawn_Qsim_YY));
for j = 1:n
    for h = 1:M
        hydro_pawn_nse{h,j} = NSE(hydro_pawn_Qsim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        hydro_pawn_pbias{h,j} = PBIAS(hydro_pawn_Qsim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        hydro_pawn_kge{h,j} = KGE(hydro_pawn_Qsim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        % Transformed objectives for minimisation
        hydro_pawn_NSEmod{h,j} = 1 - (hydro_pawn_nse{h,j}); % 1-NSE to minimise
        hydro_pawn_PBIASmod{h,j} = abs(hydro_pawn_pbias{h,j}); % abs(PBIAS) to minimise
        hydro_pawn_KGEmod{h,j} = 1 - (hydro_pawn_kge{h,j}); % 1-KGE to minimise
        % Normalise the objectives
        hydro_pawn_NSEnorm{h,j}   = (hydro_pawn_NSEmod{h,j}-min(hydro_pawn_NSEmod{h,j}))/(max(hydro_pawn_NSEmod{h,j})-min(hydro_pawn_NSEmod{h,j}));
        hydro_pawn_PBIASnorm{h,j} = (hydro_pawn_PBIASmod{h,j}-min(hydro_pawn_PBIASmod{h,j}))/(max(hydro_pawn_PBIASmod{h,j})-min(hydro_pawn_PBIASmod{h,j}));
        hydro_pawn_KGEnorm{h,j} = (hydro_pawn_KGEmod{h,j}-min(hydro_pawn_KGEmod{h,j}))/(max(hydro_pawn_KGEmod{h,j})-min(hydro_pawn_KGEmod{h,j}));
        hydro_pawn_multi{h,j} = sqrt(hydro_pawn_NSEnorm{h,j}.^2 + hydro_pawn_PBIASnorm{h,j}.^2 + hydro_pawn_KGEnorm{h,j}.^2);
        hydro_pawn_maxQ{h,j} = nanmax(hydro_pawn_Qsim_YY{h,j},[],2);
        hydro_pawn_minQ{h,j} = nanmin(hydro_pawn_Qsim_YY{h,j},[],2);
        hydro_pawn_meanQ{h,j} = nanmean(hydro_pawn_Qsim_YY{h,j},2);
    end
end

disp(' ')
%% 32. APPLY PAWN FOR THE HYDROLOGICAL MODULE
disp('APPLY PAWN FOR THE HYDROLOGICAL MODULE')

% Estimate unconditional and conditional CDFs:
% hydro_pawn_Yu = Qsim;
hydro_pawn_Yu = hydro_multi;
% hydro_pawn_YY = hydro_pawn_Qsim_YY;
hydro_pawn_YY = hydro_pawn_multi;

% [ hydro_pawn_YF, hydro_pawn_Fu, hydro_pawn_Fc  ] = pawn_cdfs(hydro_pawn_Yu,hydro_pawn_YY) ;
%
% % Plot CDFs:
% figure
% for i=1:M
%    subplot(1,M,i)
%    pawn_plot_cdf(hydro_pawn_YF, hydro_pawn_Fu, hydro_pawn_Fc(i,:),[],'y (max flow)')
% end
% 
% % Further analyze CDF of one input:
% i = 3 ;
% figure;
% pawn_plot_cdf(hydro_pawn_YF, hydro_pawn_Fu, hydro_pawn_Fc(i,:),hydro_pawn_xc{i},'y (max flow)',X_Labels{i}) % same
% % function as before but exploiting more optional input arguments
% 
% % Compute KS statistics:
% hydro_pawn_KS = pawn_ks(hydro_pawn_YF,hydro_pawn_Fu,hydro_pawn_Fc) ;
% 
% % Plot KS statistics:
% figure
% for i=1:M
%    subplot(1,M,i)
%    pawn_plot_kstest(hydro_pawn_KS(:,i),NC,NU,0.05,hydro_pawn_xc{i},X_Labels{i})
% end
% 
% % Compute PAWN index by taking a statistic of KSs (e.g. max):
% hydro_pawn_Pi = max(hydro_pawn_KS);
% 
% % Plot:
% figure 
% boxplot1(hydro_pawn_Pi,X_Labels)

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ hydro_pawn_T_m, hydro_pawn_T_lb, hydro_pawn_T_ub ] = pawn_indices(hydro_pawn_Yu,hydro_pawn_YY,stat,[],Nboot);
hydro_multi_thres = hydro_multi_sortindx(ceil(N*0.05));
% [ hydro_pawn_T_m, hydro_pawn_T_lb, hydro_pawn_T_ub ] = pawn_indices(hydro_pawn_Yu,hydro_pawn_YY,stat,[],Nboot,[],'below',hydro_multi_thres ) ;

% Plot:
figure; boxplot1(hydro_pawn_T_m,X_Labels,'pawn',hydro_pawn_T_lb,hydro_pawn_T_ub)
title('Model performance')

% % Convergence analysis:
% stat = 'max' ; % statistic to be applied to KSs
% hydro_pawn_NCb = [ NC/10 NC/2 NC ] ;
% hydro_pawn_NUb = [ NU/10 NU/2 NU ] ;
% 
% [ hydro_pawn_T_m_n, hydro_pawn_T_lb_n, hydro_pawn_T_ub_n ] = pawn_convergence( hydro_pawn_Yu, hydro_pawn_YY, stat, hydro_pawn_NUb, hydro_pawn_NCb,[],Nboot );
% hydro_pawn_NN = hydro_pawn_NUb+n*hydro_pawn_NCb ;
% figure; plot_convergence(hydro_pawn_T_m_n,hydro_pawn_NN,hydro_pawn_T_lb_n,hydro_pawn_T_ub_n,[],'no of evals',[],X_Labels)

disp(' ')
%% 33. APPLY PAWN FOR THE HYDROLOGICAL MODULE AND FOR MIN AND MAX FLOWS
disp('APPLY PAWN FOR THE HYDROLOGICAL MODULE AND FOR MIN AND MAX FLOWS')

% Estimate unconditional and conditional CDFs:
% hydro_pawn_Yu = Qsim;
hydro_pawn_maxYu = nanmax(Qsim,[],2);
hydro_pawn_minYu = nanmin(Qsim,[],2);
hydro_pawn_meanYu = nanmean(Qsim,2);
hydro_pawn_maxYY = hydro_pawn_maxQ;
hydro_pawn_minYY = hydro_pawn_minQ;
hydro_pawn_meanYY = hydro_pawn_meanQ;

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ hydro_pawn_max_T_m, hydro_pawn_max_T_lb, hydro_pawn_max_T_ub ] = pawn_indices(hydro_pawn_maxYu,hydro_pawn_maxYY,stat,[],Nboot);
[ hydro_pawn_min_T_m, hydro_pawn_min_T_lb, hydro_pawn_min_T_ub ] = pawn_indices(hydro_pawn_minYu,hydro_pawn_minYY,stat,[],Nboot);
[ hydro_pawn_mean_T_m, hydro_pawn_mean_T_lb, hydro_pawn_mean_T_ub ] = pawn_indices(hydro_pawn_meanYu,hydro_pawn_meanYY,stat,[],Nboot);

% Plot:
figure; boxplot1(hydro_pawn_max_T_m,X_Labels,'pawn',hydro_pawn_max_T_lb,hydro_pawn_max_T_ub)
title('Peak flows')
figure; boxplot1(hydro_pawn_min_T_m,X_Labels,'pawn',hydro_pawn_min_T_lb,hydro_pawn_min_T_ub)
title('Min flows')
figure; boxplot1(hydro_pawn_mean_T_m,X_Labels,'pawn',hydro_pawn_mean_T_lb,hydro_pawn_mean_T_ub)
title('Mean flows')

disp(' ')
%% 34. SETUP PAWN FOR THE RANDOM EROSION MODULE
disp('SETUP PAWN FOR THE RANDOM EROSION MODULE')

% NU = 150 ; % number of samples to estimate unconditional CDF
NC = 1000 ; % number of samples to estimate conditional CDFs
n  = 10 ; % number of conditioning points

DistrPar = cell(ME,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:ME; DistrPar{i} = [xymin(i) xymax(i) ] ; end % this is to resort to
% the format required by AAT_sampling function.

% Create input/output samples to estimate the conditional output CDFs:
[ rand_erosion_pawn_XX, rand_erosion_pawn_xc ] = pawn_sampling('lhs',ME,'unif',DistrPar,n,NC);
% Normalise the areas
for j = 1:n
    for h = 1:ME
        % counter = 0; % Just for information
        for i = 1:NC
            if sum(rand_erosion_pawn_XX{h,j}(i,1:3)) > 100
                % counter = counter + 1;
                rand_erosion_pawn_XX{h,j}(i,1:3) = 100*rand_erosion_pawn_XX{h,j}(i,1:3)/sum(rand_erosion_pawn_XX{h,j}(i,1:3));
            end
        end
        % Adjust the erosion parameter values: at = (at-as)+as
        rand_erosion_pawn_XX{h,j}(:,M+1) = sum(rand_erosion_pawn_XX{h,j}(:,[M+1 M+3]),2);
        rand_erosion_pawn_XX{h,j}(:,M+2) = sum(rand_erosion_pawn_XX{h,j}(:,[M+2 M+4]),2);
    end
end
rand_erosion_pawn_C2sim_YY = pawn_model_evaluation(myfun2,rand_erosion_pawn_XX,H,rain,evap) ;
% Only consider those values when sediment observations are present
for j = 1:n
    for h = 1:ME
        rand_erosion_pawn_C2sim_YY{h,j}(:,~(sedm>0)) = 0;
    end
end

disp(' ')
%% 35. EVALUATE PAWN STATISTICS FOR THE RANDOM EROSION MODULE
disp('EVALUATE PAWN STATISTICS FOR THE RANDOM EROSION MODULE')

% Statistics
rand_erosion_pawn_nse = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_pbias = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_kge = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_NSEmod = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_PBIASmod = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_KGEmod = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_NSEnorm = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_PBIASnorm = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_KGEnorm = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_multi = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_maxC = cell(size(rand_erosion_pawn_C2sim_YY));
rand_erosion_pawn_meanC = cell(size(rand_erosion_pawn_C2sim_YY));
for j = 1:n
    for h = 1:ME
        rand_erosion_pawn_nse{h,j} = NSE(rand_erosion_pawn_C2sim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        rand_erosion_pawn_pbias{h,j} = PBIAS(rand_erosion_pawn_C2sim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        rand_erosion_pawn_kge{h,j} = KGE(rand_erosion_pawn_C2sim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        % Transformed objectives for minimisation
        rand_erosion_pawn_NSEmod{h,j} = 1 - (rand_erosion_pawn_nse{h,j}); % 1-NSE to minimise
        rand_erosion_pawn_PBIASmod{h,j} = abs(rand_erosion_pawn_pbias{h,j}); % abs(PBIAS) to minimise
        rand_erosion_pawn_KGEmod{h,j} = 1 - (rand_erosion_pawn_kge{h,j}); % 1-KGE to minimise
        % Normalise the objectives
        rand_erosion_pawn_NSEnorm{h,j}   = (rand_erosion_pawn_NSEmod{h,j}-min(rand_erosion_pawn_NSEmod{h,j}))/(max(rand_erosion_pawn_NSEmod{h,j})-min(rand_erosion_pawn_NSEmod{h,j}));
        rand_erosion_pawn_PBIASnorm{h,j} = (rand_erosion_pawn_PBIASmod{h,j}-min(rand_erosion_pawn_PBIASmod{h,j}))/(max(rand_erosion_pawn_PBIASmod{h,j})-min(rand_erosion_pawn_PBIASmod{h,j}));
        rand_erosion_pawn_KGEnorm{h,j} = (rand_erosion_pawn_KGEmod{h,j}-min(rand_erosion_pawn_KGEmod{h,j}))/(max(rand_erosion_pawn_KGEmod{h,j})-min(rand_erosion_pawn_KGEmod{h,j}));
        rand_erosion_pawn_multi{h,j} = sqrt(rand_erosion_pawn_NSEnorm{h,j}.^2 + rand_erosion_pawn_PBIASnorm{h,j}.^2 + rand_erosion_pawn_KGEnorm{h,j}.^2);
        rand_erosion_pawn_maxC{h,j} = nanmax(rand_erosion_pawn_C2sim_YY{h,j},[],2);
        [ii,~,v] = find(rand_erosion_pawn_C2sim_YY{h,j});
        rand_erosion_pawn_meanC{h,j} = accumarray(ii,v,[],@mean);
        % rand_erosion_pawn_meanC{h,j} = nanmean(rand_erosion_pawn_C2sim_YY{h,j},2);
        clear ii v
    end
end

disp(' ')
%% 36. APPLY PAWN FOR THE RANDOM EROSION MODULE
disp('APPLY PAWN FOR THE RANDOM EROSION MODULE')

% Estimate unconditional and conditional CDFs:
% rand_erosion_pawn_Yu = C2sim;
rand_erosion_pawn_Yu = rand_erosion_multi;
% rand_erosion_pawn_YY = rand_erosion_pawn_C2sim_YY;
rand_erosion_pawn_YY = rand_erosion_pawn_multi;

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ rand_erosion_pawn_T_m, rand_erosion_pawn_T_lb, rand_erosion_pawn_T_ub ] = pawn_indices(rand_erosion_pawn_Yu,rand_erosion_pawn_YY,stat,[],Nboot);
rand_erosion_multi_thres = rand_erosion_multi_sortindx(ceil(N*0.05));
% [ rand_erosion_pawn_T_m, rand_erosion_pawn_T_lb, rand_erosion_pawn_T_ub ] = pawn_indices(rand_erosion_pawn_Yu,rand_erosion_pawn_YY,stat,[],Nboot,[],'below',rand_erosion_multi_thres ) ;

% Plot:
figure; boxplot1(rand_erosion_pawn_T_m,XY_Labels,'pawn',rand_erosion_pawn_T_lb,rand_erosion_pawn_T_ub)
title('Model performance')

disp('')
%% 37. APPLY PAWN FOR THE RANDOM EROSION MODULE AND FOR MIN AND MAX FLOWS
disp('APPLY PAWN FOR THE RANDOM EROSION MODULE AND FOR MIN AND MAX FLOWS')

% Estimate unconditional and conditional CDFs:
% rand_erosion_pawn_Yu = C2sim;
rand_erosion_pawn_maxYu = nanmax(C2sim,[],2);
rand_erosion_pawn_meanYu = nanmean(C2sim,2);
rand_erosion_pawn_maxYY = rand_erosion_pawn_maxC;
rand_erosion_pawn_meanYY = rand_erosion_pawn_meanC;

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ rand_erosion_pawn_max_T_m, rand_erosion_pawn_max_T_lb, rand_erosion_pawn_max_T_ub ] = pawn_indices(rand_erosion_pawn_maxYu,rand_erosion_pawn_maxYY,stat,[],Nboot);
[ rand_erosion_pawn_mean_T_m, rand_erosion_pawn_mean_T_lb, rand_erosion_pawn_mean_T_ub ] = pawn_indices(rand_erosion_pawn_meanYu,rand_erosion_pawn_meanYY,stat,[],Nboot);

% Plot:
figure; boxplot1(rand_erosion_pawn_max_T_m,XY_Labels,'pawn',rand_erosion_pawn_max_T_lb,rand_erosion_pawn_max_T_ub)
title('Peak sediment')
figure; boxplot1(rand_erosion_pawn_mean_T_m,XY_Labels,'pawn',rand_erosion_pawn_mean_T_lb,rand_erosion_pawn_mean_T_ub)
title('Mean sediment')

disp(' ')
%% 38. SETUP PAWN FOR THE EROSION MODULE WITH BEHAVIOURAL HYDROLOGICAL SETS
disp('SETUP PAWN FOR THE EROSION MODULE WITH BEHAVIOURAL HYDROLOGICAL SETS')

% NU = 150 ; % number of samples to estimate unconditional CDF
NC = ceil(1000/Nbeh) ; % number of samples to estimate conditional CDFs
n  = 10 ; % number of conditioning points

DistrPar = cell(E,1);
% Create data structure for parameter ranges as required by AAT_sampling
for i=1:E; DistrPar{i} = [ymin(i) ymax(i) ] ; end % this is to resort to
% the format required by AAT_sampling function.

% Create input/output samples to estimate the conditional output CDFs:
[ erosion_pawn_XX, erosion_pawn_xc ] = pawn_sampling('lhs',E,'unif',DistrPar,n,NC);
% Normalise the areas
erosion_pawn_XY = cell(E,n);
for j = 1:n
    for h = 1:E
        % Adjust the erosion parameter values: at = (at-as)+as
        erosion_pawn_XX{h,j}(:,1) = sum(erosion_pawn_XX{h,j}(:,[1 3]),2);
        erosion_pawn_XX{h,j}(:,2) = sum(erosion_pawn_XX{h,j}(:,[2 4]),2);
        % Combine behavioural hydro parameters with erosion random parameters
        erosion_pawn_XY{h,j} = nan(Nbeh*NC,M+E);
        erosion_pawn_XY{h,j}(:,M+1:M+E) = repmat(erosion_pawn_XX{h,j},Nbeh,1);
        for i = 1:Nbeh
            erosion_pawn_XY{h,j}(NC*(i-1)+1:NC*i,1:M) = repmat(hydro_multi_Xbeh(i,:),NC,1);
        end
        erosion_pawn_XX{h,j} = erosion_pawn_XY{h,j}(:,M+1:M+E);
    end
end

erosion_pawn_Csim_YY = pawn_model_evaluation(myfun2,erosion_pawn_XY,H,rain,evap) ;
% Only consider those values when sediment observations are present
for j = 1:n
    for h = 1:E
        erosion_pawn_Csim_YY{h,j}(:,~(sedm>0)) = 0;
    end
end

disp(' ')
%% 39. EVALUATE PAWN STATISTICS FOR THE BEHAVIOURAL EROSION MODULE
disp('EVALUATE PAWN STATISTICS FOR THE BEHAVIOURAL EROSION MODULE')

% Statistics
erosion_pawn_nse = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_pbias = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_kge = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_NSEmod = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_PBIASmod = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_KGEmod = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_NSEnorm = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_PBIASnorm = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_KGEnorm = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_multi = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_maxC = cell(size(erosion_pawn_Csim_YY));
erosion_pawn_meanC = cell(size(erosion_pawn_Csim_YY));
for j = 1:n
    for h = 1:E
        erosion_pawn_nse{h,j} = NSE(erosion_pawn_Csim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        erosion_pawn_pbias{h,j} = PBIAS(erosion_pawn_Csim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        erosion_pawn_kge{h,j} = KGE(erosion_pawn_Csim_YY{h,j}(:,warmup+1:end),flow(warmup+1:end)');
        % Transformed objectives for minimisation
        erosion_pawn_NSEmod{h,j} = 1 - (erosion_pawn_nse{h,j}); % 1-NSE to minimise
        erosion_pawn_PBIASmod{h,j} = abs(erosion_pawn_pbias{h,j}); % abs(PBIAS) to minimise
        erosion_pawn_KGEmod{h,j} = 1 - (erosion_pawn_kge{h,j}); % 1-KGE to minimise
        % Normalise the objectives
        erosion_pawn_NSEnorm{h,j}   = (erosion_pawn_NSEmod{h,j}-min(erosion_pawn_NSEmod{h,j}))/(max(erosion_pawn_NSEmod{h,j})-min(erosion_pawn_NSEmod{h,j}));
        erosion_pawn_PBIASnorm{h,j} = (erosion_pawn_PBIASmod{h,j}-min(erosion_pawn_PBIASmod{h,j}))/(max(erosion_pawn_PBIASmod{h,j})-min(erosion_pawn_PBIASmod{h,j}));
        erosion_pawn_KGEnorm{h,j} = (erosion_pawn_KGEmod{h,j}-min(erosion_pawn_KGEmod{h,j}))/(max(erosion_pawn_KGEmod{h,j})-min(erosion_pawn_KGEmod{h,j}));
        erosion_pawn_multi{h,j} = sqrt(erosion_pawn_NSEnorm{h,j}.^2 + erosion_pawn_PBIASnorm{h,j}.^2 + erosion_pawn_KGEnorm{h,j}.^2);
        erosion_pawn_maxC{h,j} = nanmax(erosion_pawn_Csim_YY{h,j},[],2);
        [ii,~,v] = find(erosion_pawn_Csim_YY{h,j});
        erosion_pawn_meanC{h,j} = accumarray(ii,v,[],@mean);
        clear ii v
    end
end

disp(' ')
%% 40. APPLY PAWN FOR THE EROSION MODULE WITH BEHAVIOURAL HYDROLOGICAL SETS
disp('APPLY PAWN FOR THE EROSION MODULE WITH BEHAVIOURAL HYDROLOGICAL SETS')

% Estimate unconditional and conditional CDFs:
% erosion_pawn_Yu = C2sim;
erosion_pawn_Yu = erosion_multi;
% erosion_pawn_YY = erosion_pawn_C2sim_YY;
erosion_pawn_YY = erosion_pawn_multi;

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ erosion_pawn_T_m, erosion_pawn_T_lb, erosion_pawn_T_ub ] = pawn_indices(erosion_pawn_Yu,erosion_pawn_YY,stat,[],Nboot);
erosion_multi_thres = erosion_multi_sortindx(ceil(N3*0.05));
% [ erosion_pawn_T_m, erosion_pawn_T_lb, erosion_pawn_T_ub ] = pawn_indices(erosion_pawn_Yu,erosion_pawn_YY,stat,[],Nboot,[],'below',erosion_multi_thres ) ;

% Plot:
figure; boxplot1(erosion_pawn_T_m,Y_Labels,'pawn',erosion_pawn_T_lb,erosion_pawn_T_ub)
title('Model performance')

disp('')
%% 41. APPLY PAWN FOR THE BEHAVIOURAL EROSION MODULE AND FOR MIN AND MAX FLOWS
disp('APPLY PAWN FOR THE BEHAVIOURAL EROSION MODULE AND FOR MIN AND MAX FLOWS')

% Estimate unconditional and conditional CDFs:
% erosion_pawn_Yu = C2sim;
erosion_pawn_maxYu = max(Csim,[],2);
erosion_pawn_meanYu = min(C2sim,[],2);
erosion_pawn_maxYY = erosion_pawn_maxC;
erosion_pawn_meanYY = erosion_pawn_meanC;

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ erosion_pawn_max_T_m, erosion_pawn_max_T_lb, erosion_pawn_max_T_ub ] = pawn_indices(erosion_pawn_maxYu,erosion_pawn_maxYY,stat,[],Nboot);
[ erosion_pawn_mean_T_m, erosion_pawn_mean_T_lb, erosion_pawn_mean_T_ub ] = pawn_indices(erosion_pawn_meanYu,erosion_pawn_meanYY,stat,[],Nboot);

% Plot:
figure; boxplot1(erosion_pawn_max_T_m,Y_Labels,'pawn',erosion_pawn_max_T_lb,erosion_pawn_max_T_ub)
title('Max sediment')
figure; boxplot1(erosion_pawn_mean_T_m,Y_Labels,'pawn',erosion_pawn_mean_T_lb,erosion_pawn_mean_T_ub)
title('Mean sediment')

disp(' ')
%% 42. PLOT SENSITIVITY INDICES (MVD AND PAWN) OF BEHAVIOURAL AND RANDOM SETUPS
disp('PLOT SENSITIVITY INDICES (MVD AND PAWN) OF BEHAVIOURAL AND RANDOM SETUPS')

% Plot results mvd:
figure
subplot(2,2,1)
boxplot1([hydro_multi_mvd,erosion_mvd],XY_Labels,'mvd',[hydro_multi_mvd_lb,erosion_mvd_lb],[hydro_multi_mvd_ub,erosion_mvd_ub],10)
hold on
plot([M M]+0.5,[-0.5 1.5],'--k')
grid on
title('(a) RSA of the hydrological and erosion modules')
set(gca,'XTickLabelRotation',90)

subplot(2,2,3)
boxplot1(rand_erosion_mvd,XY_Labels,'mvd',rand_erosion_mvd_lb,rand_erosion_mvd_ub,10)
grid on
title('(b) RSA of the erosion module only')
set(gca,'XTickLabelRotation',90)

% Plot results pawn:
% figure
subplot(2,2,2)
boxplot1([hydro_pawn_T_m,erosion_pawn_T_m],XY_Labels,'pawn',[hydro_pawn_T_lb,erosion_pawn_T_lb],[hydro_pawn_T_ub,erosion_pawn_T_ub],10)
hold on
plot([M M]+0.5,[-0.5 1.5],'--k')
grid on
title('(c) PAWN SA of the hydrological and erosion modules')
set(gca,'XTickLabelRotation',90)

subplot(2,2,4)
boxplot1(rand_erosion_pawn_T_m,XY_Labels,'pawn',rand_erosion_pawn_T_lb,rand_erosion_pawn_T_ub,10)
grid on
title('(d) PAWN SA of the erosion module only')
set(gca,'XTickLabelRotation',90)

disp(' ')
%% 43. PLOT SENSITIVITY INDICES (PAWN) OF OUTPUT ANALYSIS
disp('PLOT SENSITIVITY INDICES (PAWN) OF OUTPUT ANALYSIS')

% Plot results pawn:
figure
np = 3;

subplot(3,(M+ME+2)*np,(1:M*np))
boxplot1(hydro_pawn_max_T_m,X_Labels,'pawn',hydro_pawn_max_T_lb,hydro_pawn_max_T_ub,10)
grid on
title('(a) Sensitivity in peak flow simulation')
set(gca,'XTickLabelRotation',90)

subplot(3,(M+ME+2)*np,(1:M*np)+(M+ME+2)*np)
boxplot1(hydro_pawn_mean_T_m,X_Labels,'pawn',hydro_pawn_mean_T_lb,hydro_pawn_mean_T_ub,10)
grid on
title('(b) Sensitivity in mean flow simulation')
set(gca,'XTickLabelRotation',90)

subplot(3,(M+ME+2)*np,(1:M*np)+(M+ME+2)*np*2)
title('(c) Sensitivity in minimum flow simulation')
boxplot1(hydro_pawn_min_T_m,X_Labels,'pawn',hydro_pawn_min_T_lb,hydro_pawn_min_T_ub,10)
grid on
set(gca,'XTickLabelRotation',90)

subplot(3,(M+ME+2)*np,((M+ME+2)*3*np-E*np-1:(M+ME+2)*3*np)-(M+ME+2)*np*2)
boxplot1(erosion_pawn_max_T_m,Y_Labels,'pawn',erosion_pawn_max_T_lb,erosion_pawn_max_T_ub,10)
hold on
plot([M M]+0.5,[-0.5 1.5],'--k')
grid on
t = title('(d) Sensitivity in peak sediment concentration (behavioural sets)');
set(gca,'XTickLabelRotation',90)
tpox = get(t,'Position'); tpox(1) = -2.5;
set(t,'Position',tpox)

subplot(3,(M+ME+2)*np,((M+ME+2)*3*np-ME*np+1:(M+ME+2)*3*np)-(M+ME+2)*np)
boxplot1(rand_erosion_pawn_max_T_m,XY_Labels,'pawn',rand_erosion_pawn_max_T_lb,rand_erosion_pawn_max_T_ub,10)
grid on
title('(e) Sensitivity in peak sediment concentration (random sets)')
set(gca,'XTickLabelRotation',90)%,'YTick',[],'YTickLabels',[],'ylabel',[])

subplot(3,(M+ME+2)*np,((M+ME+2)*3*np-ME*np+1:(M+ME+2)*3*np))
boxplot1(rand_erosion_pawn_mean_T_m,XY_Labels,'pawn',rand_erosion_pawn_mean_T_lb,rand_erosion_pawn_mean_T_ub,10)
grid on
title('(f) Sensitivity in mean sediment concentration')
set(gca,'XTickLabelRotation',90)%,'YTick',[],'YTickLabels',[],'ylabel',[])

clear np t tpox

disp(' ')
%% 44. COORDINATE PLOTS OF BEHAVIOURAL SETS
disp('COORDINATE PLOTS OF BEHAVIOURAL SETS')

figure(135)
% Parallel coordinate plot:
subplot(2,ME*3,1:M*3+1)
parcoor(X,X_Labels,[],hydro_multi_idxb_rsa,12)
hold on
multi_bestnorm_nse = (hydro_nse_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_pbias = (hydro_pbias_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_kge = (hydro_kge_sortset(1,:)-min(X))./(max(X)-min(X));
multi_bestnorm_multi = (hydro_multi_sortset(1,:)-min(X))./(max(X)-min(X));
plot(multi_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(multi_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

% Parallel coordinate plot:
subplot(2,ME*3,ME*3+1:ME*3*2)
parcoor(XY2,XY_Labels,[],rand_erosion_idxb_rsa,10)
hold on
rand_erosion_bestnorm_nse = (rand_erosion_nse_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_pbias = (rand_erosion_pbias_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_kge = (rand_erosion_kge_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
rand_erosion_bestnorm_multi = (rand_erosion_multi_sortset(1,:)-min(XY2))./(max(XY2)-min(XY2));
plot(rand_erosion_bestnorm_nse,'-v','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_pbias,'-^','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_kge,'-s','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')
plot(rand_erosion_bestnorm_multi,'-o','MarkerSize',6,'Color',0.9*[1 1 1],'LineWidth',1,'MarkerFaceColor',0.8*[1 1 1],'MarkerEdgeColor','k')

disp('')