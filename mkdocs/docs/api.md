# API reference

## Generate and measure

::: nongauss.surface

## Plots

Every `plot_*` function returns a matplotlib figure (at most two panels), or a Plotly figure with `interactive=True`.

::: nongauss.plots
    options:
      members: [plot, KINDS, plot_surfaces, plot_heightmap_mask, plot_distribution, plot_profile_bearing, plot_hhcf_report, plot_spectrum, plot_feasibility, plot_compare, plot_logo, save]

::: nongauss.interactive
    options:
      members: [explorer]

## Johnson distributions

::: nongauss.johnson
    options:
      members: [f_johnson_M, johnson_dist, sample_johnson, johnson_boundaries]

## Command line

```bash
nongauss --skew 1 --kurt 4 --alpha 0.7 --out out/ [--interactive]
```
