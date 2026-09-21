function res = dm_difffreq(u, v, delta_u, delta_v, varargin)
%DM_DIFFFREQ  Two-sample test comparing forecasts made at different frequencies.
%
%   RES = DM_DIFFFREQ(U, V, DELTA_U, DELTA_V) tests whether two forecasts have
%   equal expected loss when their losses are observed at different sampling
%   frequencies, following Zhang, Zu and Li, "Statistical Comparison of
%   Forecasts Made at Different Frequencies".
%
%   Inputs:
%     U        n_u-by-1 loss series of forecast 1 (the HIGHER-frequency series)
%     V        n_v-by-1 loss series of forecast 2 (the LOWER-frequency series)
%     DELTA_U  sampling interval of U in common time units (e.g. 1 = monthly)
%     DELTA_V  sampling interval of V in common time units (e.g. 3 = quarterly)
%   U and V are per-period losses (e.g. squared or absolute forecast errors);
%   a positive mean(U)-mean(V) means forecast 1 has the larger average loss
%   (forecast 2 is better). Use DELTA_U <= DELTA_V (U at least as frequent).
%
%   RES = DM_DIFFFREQ(..., 'Name', Value) options:
%     'phi'       AR coefficient used in the Carlstein block-length rule
%                 (default 0.4, as in the paper's simulations)
%     'blocklen'  set the block length directly (overrides 'phi')
%     'lags'      lags for the Newey-West long-run variance of the aligned DM
%                 benchmark (default min(floor(0.75*m^(1/3)), m))
%
%   RES is a struct:
%     dm2s,  pval_dm2s   two-sample DM statistic and two-sided p-value
%                        (standard-normal reference)  <- the paper's test
%     dmcl,  pval_dmcl   clustered-t statistic and two-sided p-value
%                        (Student-t reference, nb-1 d.f.)
%     dm,    pval_dm     aligned Diebold-Mariano benchmark (normal reference)
%     nb                 number of blocks
%     blocklen           block length used
%
%   Reject equal expected loss at level ALPHA when the relevant p-value < ALPHA.
%
%   Reference:
%     Zhang, R., Zu, Y. and Li, W. Statistical comparison of forecasts made at
%     different frequencies. Journal of Business and Economic Statistics.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, 'lib'));

    ip = inputParser;
    ip.addRequired('u',  @(x) isnumeric(x) && isvector(x));
    ip.addRequired('v',  @(x) isnumeric(x) && isvector(x));
    ip.addRequired('delta_u', @(x) isscalar(x) && x > 0);
    ip.addRequired('delta_v', @(x) isscalar(x) && x > 0);
    ip.addParameter('phi', 0.4, @(x) isscalar(x) && abs(x) < 1);
    ip.addParameter('blocklen', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
    ip.addParameter('lags', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
    ip.parse(u, v, delta_u, delta_v, varargin{:});
    phi      = ip.Results.phi;
    blocklen = ip.Results.blocklen;
    lags     = ip.Results.lags;

    u = u(:);  v = v(:);
    n_u = numel(u);  n_v = numel(v);

    % ---- aligned Diebold-Mariano benchmark: subsample U onto V's timestamps ----
    u_s_align = zeros(n_v, 1);
    iu = 1;  iv = 1;
    while (iu <= n_u && iv <= n_v)
        if ((iu-1)*delta_u < (iv-1)*delta_v)
            iu = iu + 1;
        elseif ((iu-1)*delta_u == (iv-1)*delta_v)
            u_s_align(iv) = u(iu);  iv = iv + 1;
        else
            u_s_align(iv) = u(iu-1);  iv = iv + 1;
        end
    end
    d = u_s_align - v;
    if isempty(lags)
        lags = min(floor(0.75*numel(d)^(1/3)), numel(d));
    end
    dm = test_dm(d, lags);

    % ---- block construction for the two-sample and clustered-t statistics ----
    if isempty(blocklen)
        lb_u = carlstein_b(phi, n_u);
    else
        lb_u = round(blocklen);
    end
    nb = fix(n_u / lb_u);                       % number of full blocks
    lb_u_last = n_u - lb_u*nb;                  % trailing partial block
    u_bar_subs = mean(reshape(u(1:end-lb_u_last), lb_u, []), 1)';

    % aggregate V into the matching time blocks
    v_bar_subs = zeros(nb, 1);
    counter    = zeros(nb, 1);
    kb = 1;  kv = 1;
    while (kv <= n_v && kb <= nb)
        if (kv*delta_v <= kb*delta_u*lb_u)
            v_bar_subs(kb) = v_bar_subs(kb) + v(kv);
            counter(kb)    = counter(kb) + 1;
            kv = kv + 1;
        else
            kb = kb + 1;
        end
    end
    if any(counter == 0)
        error('dm_difffreq:emptyBlock', ...
            ['A time block contains no low-frequency (V) observations. ', ...
             'Increase the block length or check delta_u/delta_v.']);
    end
    v_bar_subs = v_bar_subs ./ counter;

    cov_uvbar0 = mean((u_bar_subs - v_bar_subs).^2);   % variance under the null
    cov_uvbar1 = var(u_bar_subs - v_bar_subs);         % variance under the alternative

    dm2s = (mean(u) - mean(v)) / sqrt(cov_uvbar0) * sqrt(nb);
    dmcl =  mean(u_bar_subs - v_bar_subs) / sqrt(cov_uvbar1) * sqrt(nb);

    % ---- two-sided p-values (no Statistics Toolbox required) ----
    res.dm2s      = dm2s;
    res.pval_dm2s = erfc(abs(dm2s)/sqrt(2));                  % 2*(1-normcdf(|z|))
    res.dmcl      = dmcl;
    res.pval_dmcl = t_pvalue_twosided(dmcl, nb-1);
    res.dm        = dm;
    res.pval_dm   = erfc(abs(dm)/sqrt(2));
    res.nb        = nb;
    res.blocklen  = lb_u;
end

function p = t_pvalue_twosided(tstat, df)
% Two-sided p-value of a Student-t statistic via the regularized incomplete
% beta function (base MATLAB): P(|T|>|t|) = I_{df/(df+t^2)}(df/2, 1/2).
    p = betainc(df/(df + tstat.^2), df/2, 0.5);
end
