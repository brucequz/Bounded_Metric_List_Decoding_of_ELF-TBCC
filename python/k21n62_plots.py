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

  N = 62
  K = 21
  R = K / N

  # EbNo
  ebno_dB = np.arange(1, 6.1, 0.1)
  ebno_linear = 10**(0.1*ebno_dB)
  # EsNo
  esno_linear = ebno_linear * R
  esno_dB = 10*np.log10(esno_linear)
  # Es/sigma^2
  es_over_sigma_sqrd_linear = esno_linear * 2
  es_over_sigma_sqrd_dB = 10*np.log10(es_over_sigma_sqrd_linear)

  CRC_11101101001 = dist_spectra(crc='11101101001', 
                       hamming_dist=np.array([0,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,62]), 
                       num_cwds=np.array([1,217,1457,9207,30969,87699,191301,315859,411866,411866,315859,191301,87699,30969,9207,1457,217,1]),
                       dmin=16)
  dsub_11101101001 = dsu(CRC_11101101001, esno_linear)

  # - K21N62
  SLVD_K21N62_fer = sim_result(num_sims=np.array([134732, 506099, 3026802, 18684105, 165563627]),
                                num_errors=np.array([200, 200, 200, 200, 200]),
                                ebno_dB=np.array([3, 3.5, 4, 4.5, 5]))

  # - K26N62
  SLVD_K26N62_fer = sim_result(num_sims=np.array([113368, 486360, 3086557, 18338708, 129443731]),
                                num_errors=np.array([200, 200, 200, 200, 200]),
                                ebno_dB=np.array([3, 3.5, 4, 4.5, 5]))
  
  rcu_filename = 'k21n62_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu = mat_contents['k21n62_rcu']

  # NA
  P_na = []
  for esno in esno_linear:
    p_na = normal_approx(esno, N, R)
    P_na.append(p_na)
  # np.savez("k21n62_na.npz",
  #          ebno_dB=ebno_dB,
  #          esno_linear=esno_linear,
  #          K=K,
  #          N=N,
  #          P_na=P_na)
  

  plt.figure(figsize=(7, 5.5))
  plt.semilogy(SLVD_K21N62_fer.ebno_dB, SLVD_K21N62_fer.num_errors/SLVD_K21N62_fer.num_sims, '-o', 
             linewidth=1, 
             label=r'(62, 21), SLVD($L=1e8$)', 
             markerfacecolor='none')
  dsub, = plt.semilogy(ebno_dB, dsub_11101101001, linewidth=1.5, markevery=5, label=r'$p(x)=11101101001$')
  # np.savez("k21n62_11101101001_dsub.npz", ebno_dB=ebno_dB, esno_linear=esno_linear, K=K, N=N, dsub=dsub_11101101001)

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