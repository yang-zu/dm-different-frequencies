# Comparing forecasts made at different frequencies — MATLAB

MATLAB code implementing the two-sample forecast-comparison test proposed in:

> Zhang, R., Zu, Y. and Li, W. Statistical comparison of forecasts made at different
> frequencies. *Journal of Business and Economic Statistics* (forthcoming).
> [doi:10.1080/07350015.2026.2634825](https://doi.org/10.1080/07350015.2026.2634825)

Standard Diebold–Mariano tests require the two forecasts' losses to be observed at the
**same** frequency. This test compares the average loss of two forecasts observed at
**different** frequencies (for example a monthly forecast against a quarterly one),
using a block construction to align the series in time. This is a **ready-to-use
implementation**, not replication code for the paper's tables.

## Install

```matlab
addpath(genpath('path/to/dm-different-frequencies'));
```

Requires MATLAB (no toolboxes — the p-values use base `erfc`/`betainc`).

## Usage

```matlab
% u : losses of forecast 1, the HIGHER-frequency series (e.g. monthly)
% v : losses of forecast 2, the LOWER-frequency series  (e.g. quarterly)
% delta_u, delta_v : their sampling intervals in common time units (1 and 3 here)
res = dm_difffreq(u, v, 1, 3);

fprintf('two-sample DM = %.3f  (p = %.3f)\n', res.dm2s, res.pval_dm2s);
```

`u` and `v` are per-period losses (e.g. squared or absolute forecast errors). A positive
`mean(u)-mean(v)` means forecast 1 has the larger average loss, i.e. forecast 2 is better.

`res` is a struct:

| field | meaning |
|-------|---------|
| `dm2s`, `pval_dm2s` | the paper's two-sample DM test (standard-normal reference), two-sided p-value |
| `dmcl`, `pval_dmcl` | clustered-t test (Student-t reference, `nb-1` d.f.) |
| `dm`, `pval_dm`     | aligned Diebold–Mariano benchmark (normal reference) |
| `nb`, `blocklen`    | number of blocks and the block length used |

Options: `'phi'` (AR coefficient for the Carlstein block-length rule, default `0.4`),
`'blocklen'` (set the block length directly), `'lags'` (Newey–West lags for the aligned
benchmark), and `'tail'` — `'two'` (default), `'right'` (H1: forecast 1 worse) or
`'left'` (H1: forecast 1 better), which sets whether the reported p-values are two- or
one-sided. See `help dm_difffreq`.

See [`examples/demo.m`](examples/demo.m) for a worked example. The user-facing function
is [`dm_difffreq.m`](dm_difffreq.m); supporting routines are in [`lib/`](lib/).

## How to cite

```bibtex
@article{zhang_forecasts_different_frequencies,
  title   = {Statistical comparison of forecasts made at different frequencies},
  author  = {Zhang, Rongmao and Zu, Yang and Li, Wei},
  journal = {Journal of Business and Economic Statistics},
  year    = {2026},
  doi     = {10.1080/07350015.2026.2634825}
}
```

## Notes

Extracted from the paper's Monte Carlo code and verified to reproduce it bit-for-bit;
tested on MATLAB R2025b. In a size study (n = 200) the two-sample test rejects at close
to the nominal level. The block length defaults to the Carlstein rule with AR coefficient
`phi = 0.4`, matching the paper's simulations; pass `'phi'` or `'blocklen'` to change it.

## License

MIT — see [LICENSE](LICENSE).
