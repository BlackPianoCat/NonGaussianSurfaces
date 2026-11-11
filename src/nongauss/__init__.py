"""nongauss: Gaussian and non-Gaussian self-affine surface generation.

Modules: johnson (fit + distributions), surface (generation + stats), plots, interactive.
Credits: see module headers (Yang et al. 2014; V. Constantoudis; D. Dave; M. Pierini; S. Korsak).
"""
from .surface import SurfaceResult, generate_surface, gaussian_surface, rank_map, SAimage_fft_2
from .johnson import f_johnson_M, johnson_dist, sample_johnson
from .surface import moments, hhcf_1d, abbott_curve, radial_psd, report_text, print_report

__all__ = ['SurfaceResult', 'generate_surface', 'gaussian_surface', 'rank_map', 'SAimage_fft_2',
           'f_johnson_M', 'johnson_dist', 'sample_johnson', 'moments', 'hhcf_1d', 'abbott_curve', 'radial_psd',
           'report_text', 'print_report']
