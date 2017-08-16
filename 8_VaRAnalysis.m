%% *** Part 19 - Performing VaR analysis with portfolio ***
% Parameters for VaR backtest
p                   = 0.01;            % VaR level
InvCDF              = norminv(p);
ValueAtRiskEst      = - InvCDF * ConditionalVol;

% Finding breaches
Breach              = (PortLogRet((Window+1):end, :) < - ValueAtRiskEst);
LenBreach           = size(Breach, 1);
pHat                = sum(Breach)/size(Breach, 1);
VRtest              = pHat/p;


% Unconditional coverage test
V1                  = sum(Breach);
V0                  = size(Breach, 1) - V1;

LogLConst           = V1 * log(p) + V0 * log(1-p);
LogLUnconst         = V1 * log(pHat) + V0 * log(1-pHat);
TestStatisticUC     = -2 * (LogLConst - LogLUnconst);
SignificanceUC      = 1 - chi2cdf(TestStatisticUC, 1);

% Conditional coverage test
% Estimating the transition probabilities
Aux01               = (ones(LenBreach - 1, 1) - Breach(1:(end - 1), :)) ...
                         .* Breach(2:end, :);
p01                 = sum(Aux01) / (LenBreach - 1 ...
                         - sum(Breach(1:(end -1), :)));
                     
Aux11               = Breach(1:(end - 1), :) .* Breach(2:end, :);
p11                 = sum(Aux11) / sum(Breach(1:(end -1), :));

% Computing log-likelihood; note that the constrained version here is 
% equal to the unconstrained version above
LogLConst           = LogLUnconst;

% Running loop to compute the unconstrained log-likelihood
LogLUnconst         = Breach(1, 1) * log(pHat) + ...
                         (1 - Breach(1, 1)) * log(1 - pHat);
for i=2:size(Breach, 1)
    if     (Breach(i-1, 1) == 0 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p01);
    elseif (Breach(i-1, 1) == 0 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p01);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p11);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p11);
    end
end

TestStatisticCC     = -2 * (LogLConst - LogLUnconst);
SignificanceCC      = 1 - chi2cdf(TestStatisticCC, 1);

% Calculate VaR at 1%,5% and 10%
alpha               = [0.01,0.05,0.1]; % VaR level
InvCDF              = norminv(alpha);
VaRDCC              = NaN(size(ConditionalVol,1),3);
for i = 1:3
VaRDCC(:,i)         = - InvCDF(1,i) * ConditionalVol;
end

