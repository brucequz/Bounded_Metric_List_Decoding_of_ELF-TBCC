from dataclasses import dataclass
from bounds import dsu

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
gem12_colors = [
    '#0072BD', # Blue
    '#D95319', # Orange
    '#EDB120', # Yellow
    '#7E2F8E', # Purple
    '#77AC30', # Green
    '#4DBEEE', # Light Blue
    '#A2142F', # Dark Red
    '#003E67', # Navy Blue
    '#722C0D', # Burnt Orange/Brown
    '#7C5D10', # Olive
    '#42194B', # Dark Purple
    '#3E5A19'  # Forest Green
]
plt.rcParams['axes.prop_cycle'] = cycler(color=gem12_colors)

import scipy.io

@dataclass
class dist_spectra:
  crc: str
  hamming_dist: np.array
  num_cwds: np.array
  dmin: int

@dataclass
class sim_result:
  num_sims: np.array
  num_errors: np.array
  ebno_dB: np.array

def main():
  # - K21N62
  SLVD_K21N62_fer = sim_result(num_sims=np.array([134732, 506099, 3026802, 18684105]),
                                num_errors=np.array([200, 200, 200, 200]),
                                ebno_dB=np.array([3, 3.5, 4, 4.5]))

  # - K26N62
  SLVD_K26N62_fer = sim_result(num_sims=np.array([113368, 486360, 3086557, 18338708, 129443731]),
                                num_errors=np.array([200, 200, 200, 200, 200]),
                                ebno_dB=np.array([3, 3.5, 4, 4.5, 5]))
  

  plt.figure(figsize=(7, 5.5))
  plt.semilogy(SLVD_K21N62_fer.ebno_dB, SLVD_K21N62_fer.num_errors/SLVD_K21N62_fer.num_sims, '-o', 
             linewidth=1, 
             label=r'SLVD($L=1e8$)', 
             markerfacecolor='none')
  plt.semilogy(SLVD_K26N62_fer.ebno_dB, SLVD_K26N62_fer.num_errors/SLVD_K26N62_fer.num_sims, '-s', 
             linewidth=1, 
             label=r'SLVD($L=1e8$)', 
             markerfacecolor='none')
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([2.5, 5])
  plt.ylim([1e-8, 1e-1])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel('Frame Error Rate', fontsize=15)
  plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()


if __name__ == "__main__":
  main()