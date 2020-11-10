function [x, loss] = Robust_ETF(y, lam, k)
% This function realizes Robust Enhanced Trend Filtering (Robust_ETF)
%% Input %%%%%%%%%%
%   y      : noisy data
%   lam    : trade-off parameter
%   k      : the k order difference matrix
%% Output %%%%%%%%%%
%   x      : extracted trend signal
%   loss   : the corresponding loss function

% Reference: 'Robust Enhanced Trend Filtering with Unknown Noise',
% Signal Processing, 2020
% https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2020.6


%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%
n = length(y);
D = diff(eye(n), k, 1);
M = 5;
R = initialization(y(:)',M);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));
model.mu = zeros(M, 1);
model.Sigma = rand(M, 1);
nk = sum(R,1);
model.weight = nk(:)/size(R,1);
converged = false;
c = 1;
[x, ~] = update_x(y, y, lam, ones(length(y),1), D, c);
x_old = x;
Error = y - x(:);
%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%

t = 0;
maxiter = 1000;
tol = 1.0e-8;
%%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
[R, ~] = expectation(Error, model);
%%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%% 
while ~converged &&  t < maxiter  
    t = t+1;
    
    %%%%%%%%%%%%%%%% M Step 1 %%%%%%%%%%%%%%%%%%%
    [model] = maximizationModel(Error, R);
    %%%%%%%%%%%%%%%% M Step 1 %%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% M Step 2 %%%%%%%%%%%%%%%%%%%
    [W, x, cost] = maximizationW(model, y, x, R, lam, D, c);
    %%%%%%%%%%%%%%%% M Step 2 %%%%%%%%%%%%%%%%%%%
    

    Error = y - x(:);
    %%%%%%%%%%%%%%%% Reduce M %%%%%%%%%%%%%%%%%%%
    % model = reducek(model);
    %%%%%%%%%%%%%%%% Reduce M %%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    [R, ~] = expectation(Error, model);
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%  
    
    %%%%%%%%%%%%%%%% Update c %%%%%%%%%%%%%%%%%%%
    c = max(1.0e-8, c / 10);
    %%%%%%%%%%%%%%%% Update c %%%%%%%%%%%%%%%%%%%
    
    tmp = R * (1 ./ (model.Sigma) / 2);
    loss(t) = -cost + sum(R * (log(model.weight) - log(sqrt(2*pi*model.Sigma)))) - tmp' * Error.^2;    
    
    if t >= 2
        [~,label(:)] = max(R,[],2);
        u = unique(label);   % non-empty components
        if size(R,2) ~= size(u,2)
            R = R(:,u);   % remove empty components
        else
            converged = (abs(loss(t)-loss(t-1)) < tol*abs(loss(t))) | (sqrt(sum((x-x_old).^2)) < tol*sqrt(sum(x.^2)));
        end
    end
    
    x_old = x;

end

function model = reducek(model)
mu = model.mu;
if length(mu) > 1
    Sigma = model.Sigma;
    w = model.weight;

    k = length(w);

    relative_deviation = [];
    relation = [];
    for i = 1:k-1
        relative_deviation_tmp = abs(Sigma(i) - Sigma(i+1:k)) ./ (Sigma(i) + Sigma(i+1:k));
        tmp = i+1:k;
        relation_tmp = [repmat(i, k-i, 1), tmp(:)];
        relative_deviation = [relative_deviation; relative_deviation_tmp];
        relation = [relation; relation_tmp];
    end
    [value, index] = min(relative_deviation);
    if value < 0.1
        reduce_relation = relation(index,:);
        new_weight = w(reduce_relation(1)) + w(reduce_relation(2));
        new_Sigma = (Sigma(reduce_relation(1)) + Sigma(reduce_relation(2))) / 2;
        new_index = setdiff(1:k, reduce_relation);
        Sigma = [Sigma(new_index); new_Sigma];
        w = [w(new_index); new_weight];
        mu = mu(1:k-1);
        model.mu = mu;
        model.Sigma = Sigma;
        model.weight = w;
    end
end


function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
    while k ~= length(u)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
        [u,~,label] = unique(label);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d  %initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end


function [R, llh] = expectation(y, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(y,1);
k = size(mu,1);
logRho = zeros(n,k);

for i = 1:k
    logRho(:,i) = loggausspdf(y(:)', mu(i), Sigma(i));
end
logRho = bsxfun(@plus, logRho, log(w)');
T = logsumexp(logRho, 2);
llh = sum(T) / n; % loglikelihood
logR = bsxfun(@minus, logRho, T);
R = exp(logR);


function g = loggausspdf(y, mu, Sigma)
y = y - mu;
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\y;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = log(2*pi)+2*sum(log(diag(U)));   % normalization constant
g = -(c+q)/2;


function [model] = maximizationModel(y, R)

k = size(R, 2);

nk = sum(R, 1);
mu = zeros(k, 1);%fix mu to zero
w = nk / size(R, 1);
Sigma = zeros(k, 1);
for i = 1:k
    Sigma(i) = R(:,i)' * (y - mu(i)).^2 / nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end

model.mu = mu(:);
model.Sigma = Sigma(:);
model.weight = w(:);



function [W, x, cost] = maximizationW(model, y, x, R, lam, D, c)
Sigma = model.Sigma;
W = R * (1./Sigma);
[x, cost] = update_x(y, x, lam, W, D, c);




