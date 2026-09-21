function lrv=lrvarnw(data,lags)
% Long-run covariance estimation using Newey-West (Bartlett)
% weights 
% data should be already demeaned 
%  
% USAGE:
%   lrv = lrvarnw(DATA)
%
% INPUTS:
%   data   - T by K vector of dependent data
%   lags   - Non-negative integer containing the lag length to use.  If empty or not included,
%              NLAG=min(floor(1.2*T^(1/3)),T) is used 
%
% OUTPUTS:
%   lrv      - A K by K covariance matrix estimated using Newey-West (Bartlett) weights
%   
 

T=size(data,1);
 
% NW weights
w=(lags+1-(0:lags))./(lags+1);
% Start the covariance
lrv=data'*data/T;
for i=1:lags
    Gammai=(data((i+1):T,:)'*data(1:T-i,:))/T;
    GplusGprime=Gammai+Gammai';
    lrv=lrv+w(i+1)*GplusGprime;
end


