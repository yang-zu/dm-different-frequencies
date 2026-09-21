function dmstat = test_dm(d,lags)
% Diebold Mariano test statistic
% long run variance estimator is Newey-West with Barlett 
% kernel
% input:
% d: forecast loss differential
% lags: lags used in the long-run variance estimator

n = length(d);
dmstat = sqrt(n)*mean(d)/sqrt(lrvarnw(d,lags));

end