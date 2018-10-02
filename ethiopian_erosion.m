function [C,Q,FLUXES]=ethiopian_erosion(param,H,FLUXES,hydro,rain,evap,varargin)
%
% This function runs a rainfall-runoff model developed for the Ethiopian
% highlands by Collick et al. (2009), Steenhuis et al. (2009), Tesemma et
% al. (2010), and an erosion model by Tilahun et al. (2013a, 2013b, 2015).
%
% [C,Q,OP:FLUXES] = ethiopian_erosion(param,H,FLUXES,OP:hydro,rain,evap)
% 
% Input:
%  param = vector of erosion model parameters                - vector (1,4)
%            1. at1  = Sediment transport limit for A1 [(g/L)/(mm/d)^-0.4]
%            2. at2  = Sediment transport limit for A2 [(g/L)/(mm/d)^-0.4]
%            3. as1  = Sediment source limit for A1 [(g/L)/(mm/d)^-0.4]
%            4. as2  = Sediment source limit for A2 [(g/L)/(mm/d)^-0.4]
%  H = time series of the fraction of the runoff-producing area with active
%       rill formation.
%       H ranges from 1 when rill erosion is active, to 0 when it stops.
%  FLUXES = time series of runoff, baseflow, and interflow  - vector (T,4+)
%       colum 1: RF1 = Runoff component of A1 [mm/Dt]
%       colum 2: RF2 = Runoff component of A2 [mm/Dt]
%       colum 3: BF  = Baseflow component of A3 [mm/Dt]
%       colum 4: IF  = Interflow component of A3 [mm/Dt]
%
% Optional input: If FLUXES is unknown.
%  hydro = vector of hydrological model parameters           - vector (1,9)
%            1. A1  = Percentage of saturated area [%]
%            2. A2  = Percentage of degraded area [%]
%            3. A3  = Percentage of infiltration area [%]
%                     A1+A2+A3 = 100 (%)
%            4. Sm1 = Maximum water storage A1 [mm]
%            5. Sm2 = Maximum water storage A1 > A2 [mm]
%            6. Sm3 = Field capacity of A3, generally > A1 > A2 [mm]
%            7. BSm = Maximum baseflow storage of the lineal reservoir [mm]
%            8. tm  = Half-life time of the baseflow linear reservoir [day]
%            9. tau = Duration of the interflow zero-order reservoir [day]
%  rain = time series of rainfall                            - vector (T,1)
%  evap = time series of potential evaporation               - vector (T,1)
%
%
% Output:
%  C = time series of simulated sediment                     - vector (T,1)
%  Q = time series of total flow                             - vector (T,1)
%
% Optional output: If the hydrological model was run.
%  FLUXES = time series of simulated fluxes (all in mm/Dt)   - matrix (T,9)
%
% See ethiopian_sim for references about the hydrological model.
%
% References:
%
% Ciesiolka, C. A. A., K. J. Coughlan, C. W. Rose, M. C. Escalante, G. M.
% Hashim, E. P. Paningbatan Jr., and S. Sombatpanit. 1995. Methodology for
% a multi-country study of soil erosion management. Soil Tech. 8(3): 179-192.
%
% Tilahun SA, Mukundan R, Demisse BA, Engda TA, Guzman CD, Tarakegn BC,
% Easton ZM, Collick AS, Zegeye AD, Schneiderman EM, Parlange JY, Steenhuis
% TS. 2013b. A saturation excess erosion model. Transactions of the
% American Society of Agricultural and Biological Engineers 56: 681–695.
% DOI: 10.13031/2013.42675
%
% Yu, B., C. W. Rose, B. C. Ciesiolka, K. J. Coughlan, and B. Fentie. 1997.
% Toward a framework for runoff and soil loss prediction using GUEST
% technology. Australian J. Soil Res. 35(5): 1191-1212.
%
% This function has been edited by Boris Ochoa-Tocachi in 2017 based on the
% SAFE Toolbox by F. Pianosi, F. Sarrazin and T. Wagener (2015).

% -----------------------
% Initialize variables:
% -----------------------

if size(FLUXES,2) < 4
    % disp('Running the hydrological model to generate flux data for erosion')
    M = 9; % number of hydrological model parameters
    hydro = hydro(:);
    if ~isnumeric(hydro); error('input argument ''hydro'' must be numeric'); end
    if length(hydro)~=M; error('input argument ''hydro'' must have %d components',M); end
    % Run the hydrological model:
    [~,~,FLUXES] = ethiopian_same(hydro,rain,evap);
end

M = 4; % number of erosion model parameters
param = param(:);
if ~isnumeric(param); error('input argument ''param'' must be numeric'); end
if length(param)~=M; error('input argument ''param'' must have %d components',M); end

% --------------------------
% Recover model parameters:
% --------------------------
A1  = hydro(1)/100;
A2  = hydro(2)/100;
% A3  = hydro(3)/100;
at1 = param(1); % Sediment transport limit for A1 [(g/L)/(mm/d)^-0.4]
at2 = param(2); % Sediment transport limit for A2 [(g/L)/(mm/d)^-0.4]
as1 = param(3); % Sediment source limit for A1 [(g/L)/(mm/d)^-0.4]
as2 = param(4); % Sediment source limit for A2 [(g/L)/(mm/d)^-0.4]
RF1 = FLUXES(:,1); % Runoff component of A1 [mm/day]
RF2 = FLUXES(:,2); % Runoff component of A2 [mm]
BF  = FLUXES(:,3); % Baseflow component of A3 [mm]
IF  = FLUXES(:,4); % Interflow component of A3 [mm]

n = 0.4; % Exponent value

% --------------------------
% Sediment concentration:
% --------------------------

FLUXES = [RF1 RF2 BF IF]; % Fluxes output vector
Y = A1*((RF1/A1).^(1+n)).*(as1+H*(at1-as1)) + A2*((RF2/A2).^(1+n)).*(as2+H*(at2-as2));% Sediment yield
Q = RF1 + RF2 + (BF+IF); % Total flow
C =  Y./Q;