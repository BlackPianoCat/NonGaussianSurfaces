# Johnson distributions

Module `nongauss.johnson`: turn (mean, rms, skewness, kurtosis) into a Johnson distribution and sample from it.
The low-level `sub_*` routines are the translated fitting algorithm; most users only need `f_johnson_M`,
`johnson_dist` and `sample_johnson`.

::: nongauss.johnson
    options:
      members: [f_johnson_M, johnson_dist, sample_johnson, johnson_boundaries]
