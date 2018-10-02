function [C,Q,FLUXES] = ethiopian_ped(param,H,rain,evap)
%
% This function runs the complete PED model developed for the Ethiopian
% highlands by Collick et al. (2009), Steenhuis et al. (2009), Tesemma et
% al. (2010), and an erosion model by Tilahun et al. (2013a, 2013b, 2015).
%
% [C,Q,FLUXES] = ethiopian_ped(param,H,rain,evap)
% 
% Input:
%  param = vector of PED model parameters               - vector (1,13)
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
%           10. at1  = Sediment transport limit for A1 [(g/L)/(mm/d)^-0.4]
%           11. at2  = Sediment transport limit for A2 [(g/L)/(mm/d)^-0.4]
%           12. as1  = Sediment source limit for A1 [(g/L)/(mm/d)^-0.4]
%           13. as2  = Sediment source limit for A2 [(g/L)/(mm/d)^-0.4]
%  H = time series of the fraction of the runoff-producing area with active
%       rill formation.
%       H ranges from 1 when rill erosion is active, to 0 when it stops.
%  rain = time series of rainfall                            - vector (T,1)
%  evap = time series of potential evaporation               - vector (T,1)
%
%
% Output:
%  C = time series of simulated sediment                     - vector (T,1)
%  Q = time series of total flow                             - vector (T,1)
%  FLUXES = time series of simulated fluxes (all in mm/Dt)   - matrix (T,8)
%
% See ethiopian_same for references about the hydrological model.
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

M = 13; % number of PED Model parameters
param = param(:);
if ~isnumeric(param); error('input argument ''param'' must be numeric'); end
if length(param)~=M; error('input argument ''param'' must have %d components',M); end

% --------------------------
% Recover model parameters:
% --------------------------
hydro = param(1:9);
A1  = hydro(1)/100;
A2  = hydro(2)/100;
at1 = param(10); % Sediment transport limit for A1 [(g/L)/(mm/d)^-0.4]
at2 = param(11); % Sediment transport limit for A2 [(g/L)/(mm/d)^-0.4]
as1 = param(12); % Sediment source limit for A1 [(g/L)/(mm/d)^-0.4]
as2 = param(13); % Sediment source limit for A2 [(g/L)/(mm/d)^-0.4]
n = 0.4; % Exponent value

% --------------------------
% Run the hydrological model:
% --------------------------

[~,~,FLUXES] = ethiopian_same(hydro,rain,evap);
RF1 = FLUXES(:,1); % Runoff component of A1 [mm/day]
RF2 = FLUXES(:,2); % Runoff component of A2 [mm]
BF  = FLUXES(:,3); % Baseflow component of A3 [mm]
IF  = FLUXES(:,4); % Interflow component of A3 [mm]

% --------------------------
% Sediment concentration:
% --------------------------

FLUXES = [RF1 RF2 BF IF]; % Fluxes output vector
Y = A1*((RF1/A1).^(1+n)).*(as1+H*(at1-as1)) + A2*((RF2/A2).^(1+n)).*(as2+H*(at2-as2));% Sediment yield
Q = RF1 + RF2 + (BF+IF); % Total flow
C =  Y./Q;