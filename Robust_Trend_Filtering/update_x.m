function [x, cost] = update_x(y, x, lam, W, D, c)
% This function realizes update the estimated x
%% Input %%%%%%%%%%
%   y      : noisy data
%   x      : esimated x at the previous iteration
%   lam    : trade-off parameter
%   W      : the weighted matrix
%   D      : the difference matrix
%   c      : the smoothed parameter at zero
%% Output %%%%%%%%%%
%   x      : estimated x at this itetation

% Reference: 'Robust Enhanced Trend Filtering with Unknown Noise',
% Signal Processing, 2020
% https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2020.6
W_1 = 1 ./ W;
W_1 = diag(W_1);
y = y(:);                                          

pen = 'Lp';
[phi, wfun] = penalty_funs(pen, lam, c);


DWDT = D * W_1 * D';
Dy = D*y;
Dx = D*x;
W_2 = wfun(Dx);
W_2 = diag(W_2);
F = W_2 + DWDT;                % F : Sparse banded matrix    
x = y - W_1*D'*(F\Dy);         % Solve banded linear system
x = x(:);
cost = sum(phi(D*x));          % Save cost function history


function [phi, wfun] = penalty_funs(pen, lam, c)
% [phi, dphi, wfun, ds] = penalty_funs(pen)
%
% Penalty and associated functions
%
% Input
%   pen : 'L1', 'mcp', 'Lp'
%
% Ouput
%   phi : penalty function
%   wfun : x/dphi(x)


switch pen
    case 'L1'        
        phi = @(x) sqrt(x.^2+c);
        wfun = @(x) sqrt(x.^2+c);
        
    case 'mcp'
        gamma = 1.4;
        phi = @(x) (lam .* (sqrt(x.^2+c) - (x.^2+c) ./ (2*gamma*lam))) .* (sqrt(x.^2+c) <= gamma*lam) ...
            + (gamma*lam.^2/2) .* (sqrt(x.^2+c) > gamma*lam);
        wfun = @(x) sqrt(x.^2+c) ./ (lam .* (1 - sqrt(x.^2+c) ./ (gamma*lam))) .* (sqrt(x.^2+c) <= gamma*lam) ...
            + sqrt(x.^2+c) ./ (1.0e-6) .* (sqrt(x.^2+c) > gamma*lam);

    case 'Lp'
        p = 0.5;
        phi = @(x) lam * sum((x.^2+c).^(p/2));
        wfun = @(x) ((x.^2 + c) .^ (1-p/2)) / (lam * p);   % Note: divide by zero is infinity 
end


