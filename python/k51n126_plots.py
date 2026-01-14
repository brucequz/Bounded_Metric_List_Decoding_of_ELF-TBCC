from dataclasses import dataclass
from bounds import dsu, normal_approx

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

  N = 126
  K = 51
  R = K / N

  # EbNo
  ebno_dB = np.arange(1, 5.1, 0.1)
  ebno_linear = 10**(0.1*ebno_dB)
  # EsNo
  esno_linear = ebno_linear * R
  esno_dB = 10*np.log10(esno_linear)
  # Es/sigma^2
  es_over_sigma_sqrd_linear = esno_linear * 2
  es_over_sigma_sqrd_dB = 10*np.log10(es_over_sigma_sqrd_linear)
  
  rcu_filename = 'k51n126_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu = mat_contents['k51n126_rcu']

  # NA
  P_na = []
  for esno in esno_linear:
    p_na = normal_approx(esno, N, R)
    P_na.append(p_na)
  

  plt.figure(figsize=(7, 5.5))
  # - random coding bound
  rcu, = plt.semilogy(rcu[:,1], rcu[:,3], color='k', linestyle='-.', linewidth=1.5, label='Random Coding Union Bound')

  # - normal approximation
  na, = plt.semilogy(ebno_dB, P_na, linestyle='--', linewidth=1.5, label='Normal Approximation')

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([1, 5])
  plt.ylim([1e-8, 1e-1])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel('Frame Error Rate', fontsize=15)
  plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()


if __name__ == "__main__":
  main()