function [Q_sim,STATES,FLUXES] = ethiopian_same(param,rain,evap)
%
% This function simulates a rainfall-runoff model developed for the
% Ethiopian highlands by Collick et al. (2009), Steenhuis et al. (2009),
% and Tesemma et al. (2010). This has been applied by Engda et al. (2012),
% Tilahun et al. (2013a, 2013b, 2015), Guzman et al. (2016, 2017), Zimale
% et al. (2017).
%
% Usage:
%
% [Q_sim,STATES,FLUXES] = ethiopian_sim(param,prec,evap)
%
% Input:
%  param = vector of model parameters - vector (1,9):
%            1. A1  = Percentage of saturated area [%]
%            2. A2  = Percentage of degraded area [%]
%            3. A3  = Percentage of infiltration area [%]
%                     A1+A2+A3 = 100 (%) (see line 82)
%            4. Sm1 = Maximum water storage A1 [mm]
%            5. Sm2 = Maximum water storage A2 < A1 [mm]
%            6. Sm3 = Field capacity of A3, generally > A2 > A1 [mm]
%            7. BSm = Maximum baseflow storage of the lineal reservoir [mm]
%            8. tm  = half-life time of the baseflow linear reservoir [day]
%            9. tau = duration of the interflow zero-order reservoir [day]
%  rain = time series of rainfall              - vector (T,1)
%  evap = time series of potential evaporation - vector (T,1)
%
% Output:
%  Q_sim = time series of simulated flow       - vector (T,1)
% STATES = time series of simulated states     - matrix (T,5)
% FLUXES = time series of simulated fluxes     - matrix (T,9)
%          (see comments in the code for the definition of state and flux
%          variables)
% FLUXES needed for the erosion model 'ethiopian_erosion.m':
%       colum 1: RF1 = Runoff component of A1 [mm/Dt]
%       colum 2: RF2 = Runoff component of A2 [mm/Dt]
%       colum 3: BF  = Baseflow component of A3 [mm/Dt]
%       colum 4: IF  = Interflow component of A3 [mm/Dt]
% 
% Recommended parameter values:  Smax should fall in [0,10] (mm)
%            1. A1  should fall in [0,100] (%)
%            2. A2  should fall in [0,100] (%)
%            3. A3  should fall in [0,100] (%)
%               A1+A2+A3 should sum 100 (%)
%            4. Sm1 should fall in [0,20] (mm)
%            5. Sm2 should fall in [27,300] (mm)
%            6. Sm3 should fall in [0,500] (mm)
%            7. BSm should fall in [0,500] (mm)
%            8. tm  should fall in [1,365] (day)
%            9. tau should fall in [1,365] (day)
% 
% References:
%
% Collick, A. S., Z. M. Easton, T. Ashagrie, B. Biruk, S. A. Tilahun, E.
% Adgo, S. B. Awulachew, G. Zeleke, and T. S. Steenhuis. 2009. A simple
% semidistributed water balance model for the Ethiopian highlands. Hydrol.
% Proc. 23(26): 3718-3727.
%
% Steenhuis, T. S., A. S. Collick, Z. M. Easton, E. S. Leggesse, H. K.
% Bayabil, E. D. White, S. B. Awulachew, E. Adgo, and A. A. Ahmed. 2009.
% Predicting discharge and erosion for the Abay (Blue Nile) with a simple
% model. Hydrol. Proc. 23(26): 3728-3737.
%
% Tesemma, Z. K., Y. A. Mohamed, and T. S. Steenhuis. 2010. Trends in
% rainfall and runoff in the Blue Nile basin: 1964-2003. Hydrol. Proc.
% 24(25): 3747-3758.
%
% Tilahun, S. A., C. D. Guzman, A. D. Zegeye, T. A. Engda, A. S. Collick,
% A. Rimmer, and T. S. Steenhuis. 2012. An efficient semi-distributed
% hillslope erosion model for the sub-humid Ethiopian highlands. Hydrol.
% Earth Syst. Sci. 9(2): 2121-2155.
%
% Tilahun SA, Mukundan R, Demisse BA, Engda TA, Guzman CD, Tarakegn BC,
% Easton ZM, Collick AS, Zegeye AD, Schneiderman EM, Parlange JY, Steenhuis
% TS. 2013b. A saturation excess erosion model. Transactions of the
% American Society of Agricultural and Biological Engineers 56: 681–695.
% DOI: 10.13031/2013.42675
%
% Tilahun, S.A., Guzman, C.D., Zegeye, A.D., Dagnew, D.C., Collick, A.S.,
% Yitaferu, B., Steenhuis, T.S., 2015. Distributed discharge and sediment
% concentration predictions in the sub-humid Ethiopian highlands: the Debre
% Mawi watershed. Hydrol. Process. 29, 1817–1828. DOI: 10.1002/hyp.10298. 
%
% Guzman CD, Zimale FA, Tebebu TY, et al., 2017. Modeling discharge and
% sediment concentrations after landscape interventions in a humid monsoon
% climate: The Anjeni watershed in the highlands of Ethiopia. Hydrological
% Processes, 31:1239–1257. DOI:10.1002/hyp.11092
%
%
% This function has been edited by Boris Ochoa-Tocachi in 2017 based on the
% SAFE Toolbox by F. Pianosi, F. Sarrazin and T. Wagener (2015).

% --------------------------
% Recover model parameters:
% --------------------------
if sum(param(1:3)) > 100
    A1    = 100*param(1)/sum(param(1:3)); % Fraction of saturated area [-]
    A2    = 100*param(2)/sum(param(1:3)); % Fraction of degraded area [-]
    A3    = 100*param(3)/sum(param(1:3)); % Fraction of infiltration area [-]
else
    A1    = param(1); % Fraction of saturated area [-]
    A2    = param(2); % Fraction of degraded area [-]
    A3    = param(3); % Fraction of infiltration area [-]
end
Sm1   = max(eps,param(4)) ; % Maximum Soil Moisture 1 (cannot be zero!)
Sm2   = max(eps,param(5)) ; % Maximum Soil Moisture 2 (cannot be zero!)
Sm3   = max(eps,param(6)) ; % Field Capacity for A3 (cannot be zero!)
BSm   = max(eps,param(7)) ; % Maximum baseflow storage (cannot be zero!)
alpha = log(2)/param(8) ; % Baseflow reservoir coefficient [1/day]
tau   = param(9) ; % Interflow duration [day]

% -----------------------
% Initialize variables:
% -----------------------

N = length(rain); % Number of time samples

base = zeros(N+1,1); % Baseflow reservoir moisture [mm]
perc = zeros(N,1); % Percolation in area 3 [mm/Dt]
BF = zeros(N,1); % Baseflow [mm/Dt]
IF = zeros(N+ceil(tau),1); % Interflow [mm/Dt]

Sini1 = Sm1/2; % Initial Soil Moisture 1 [mm]
Sini2 = Sm2/2; % Initial Soil Moisture 2 [mm]
Sini3 = Sm3/2; % Initial Soil Moisture 3 [mm]
base(1) = BSm/2; % Initial Baseflow Reservoir Level [mm]

% --------------------------
% Soil Moisture Dynamics:
% --------------------------

[runoff1,soil1,Ea1] = soilmoisture(rain,evap,Sm1,Sini1);
[runoff2,soil2,Ea2] = soilmoisture(rain,evap,Sm2,Sini2);
[rech,soil3,Ea3] = soilmoisture(rain,evap,Sm3,Sini3);

Ea = (A1*Ea1 + A2*Ea2 + A3*Ea3)/100; % Total actual evapotranspiration
RF1 = (A1*runoff1)/100; % Runoff produced by area 1
RF2 = (A2*runoff2)/100; % Runoff produced by area 2
RF = RF1 + RF2; % Total runoff is produced by areas 1 and 2
rech = (A3*rech)/100; % Recharge is only produced beneath area 3

% -------------------------
% Groundwater Dynamics:
% -------------------------

for t=1:N
    
    % Baseflow storage
    base(t+1) = rech(t)+base(t);
    % Add the recharge to the baseflow storage
    if base(t+1) > BSm
        perc(t)  = base(t+1) - BSm; % Compute the value of the percolation
        base(t+1) = BSm; % Limit the baseflow storage to the maximum capacity
    end
    % Baseflow dynamics
    BF(t) = base(t+1)*(1-exp(-alpha));
    base(t+1) = base(t+1) - BF(t);
    
    % Interflow dynamics
    if perc(t) ~= 0
        for i = 0:tau-1
            IF(t+i) = IF(t+i) + 2*perc(t)*(1/tau - (i+0.5)/tau^2);
            % (i+0.5) is the form used in the Excel PED model
        end
    end
end
IF(N+1:end) = [];

% -------------------------
% Model outcomes:
% -------------------------

GF = BF + IF; % Subsurface flow
Q_sim = RF + GF; % Total flow

STATES=[soil1,soil2,soil3,base];
FLUXES=[RF1,RF2,BF,IF,Ea,RF,GF,rech,perc];


function [outflow,soil,Ea] = soilmoisture(rain,evap,smax,sini)

% -----------------------
% Initialize variables:
% ----------------------

N = length(rain); % Number of time samples

soil = zeros(N+1,1); % Soil Moisture [mm]
Ea = zeros(N,1); % Actual Evapotranspiration [mm/Dt]
outflow = zeros(N,1); % Runoff in area 1 [mm/Dt]

soil(1) = sini; % Initial Soil Moisture [mm]

% --------------------------
% Soil Moisture Dynamics:
% --------------------------

for t=1:N
    
    if rain(t) < evap(t)
        soil(t+1) = soil(t)*exp((rain(t)-evap(t))/smax);
        % Add the water balance to the soil moisture content
        Ea(t) = soil(t)-soil(t+1); % Compute the actual evapotranspiration
    else
        soil(t+1) = rain(t)-evap(t)+soil(t);
        % Add the water balance to the soil moisture content
        Ea(t) = evap(t); % Compute the actual evapotranspiration
    end
    
    if soil(t+1) > smax
        outflow(t)  = soil(t+1) - smax; % Compute the value of the outflow
        soil(t+1) = smax; % Limit the soil moisture to the maximum capacity
    end
    
end