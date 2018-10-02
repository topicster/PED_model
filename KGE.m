function [kge,r,a,b] = KGE(y_sim,y_obs,scale,varargin)
%
% Computes the Kling-Gupta Efficiency (KGE) coefficient
%
% [kge,r,a,b] = KGE(Y_sim,y_obs)
% 
% y_sim = time series of modelled variable     - matrix (N,T)
%         (N>1 different time series can be evaluated at once)
% y_obs = time series of observed variable     - vector (1,T)
% scale = scaling factors for r, a, and b      - vector (1,3)
%
% kge   = vector of KGE coefficients           - vector (N,1)
% r     = linear correlation coefficient between y_sim and y_obs - vector (N,1)
% a     = relative variability as the ratio of std of y_sim and y_obs - vector (N,1)
% b     = bias as the ratio of means between y_sim and y_obs - vector (N,1)
%
% References:
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., and Martinez, G. F. (2009):
% Decomposition of the mean squared error and NSE performance criteria:
% Implications for improving hydrological modelling, J. Hydrol., 377,
% 80–91, doi:10.1016/j.jhydrol.2009.08.003

% This function has been modified by Boris Ochoa-Tocachi, based on a script
% part of the SAFE Toolbox by F. Pianosi, F. Sarrazin and T. Wagener at
% Bristol University (2015).
% boris.ochoa13@imperial.ac.uk
%
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

% Check the input variables
if nargin < 3
    scale = [1 1 1];
elseif length(scale) ~= 3
    error('input scale must have 3 parameters.')
end
if nargin < 2
    error('Not enough input arguments.')
end

[N,T] = size(y_sim) ;
[M,D] = size(y_obs) ;

if T~=D
    error('input y_sim and y_obs must have the same number of columns.')
end
if M>1
    error('input y_obs must be a row vector.')
end

% Bias
muO = mean(y_obs);
muS = mean(y_sim,2);
b = muS / muO;

% Variability
stdO = std(y_obs,1);
stdS = std(y_sim,1,2);
a = stdS / stdO;

% Correlation
ErrO = (y_obs - muO) / stdO;
ErrS = (y_sim - muS) ./ stdS;
Err  = repmat(ErrO,N,1) .* ErrS;

r = (1/T)*sum(Err,2);

kge = 1 - sqrt(scale(1)*(r-1).^2 + scale(2)*(a-1).^2 + scale(3)*(b-1).^2);