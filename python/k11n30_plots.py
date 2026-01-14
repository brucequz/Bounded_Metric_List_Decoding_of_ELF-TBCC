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
  hamming_distance = np.array([0,6,8,10,12,14,16,18,20,22,24,30])
  

  ELF_10001 = dist_spectra(crc='10001', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,4,20,97,358,541,571,299,138,19,0,0]),
                       dmin=6)
  
  CRC_10011 = dist_spectra(crc='10011', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,0,30,108,300,585,585,300,108,30,0,1]),
                       dmin=8)
  
  ELF_10101 = dist_spectra(crc='10101', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,0,20,103,360,561,511,361,100,31,0,0]),
                       dmin=8)
  
  ELF_10111 = dist_spectra(crc='10111', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,3,23,107,335,567,512,365,121,14,0,0]),
                       dmin=6)
  
  CRC_11001 = dist_spectra(crc='11001', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,45,45,258,270,405,405,270,258,45,45,1]),
                       dmin=6)
  
  ELF_11011 = dist_spectra(crc='11011', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,0,21,100,367,544,534,348,101,32,0,0]),
                       dmin=8)
  
  ELF_11101 = dist_spectra(crc='11101', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,4,13,107,343,589,512,333,121,23,2,0]),
                       dmin=8)
  
  CRC_11111 = dist_spectra(crc='11111', 
                       hamming_dist=hamming_distance, 
                       num_cwds=np.array([1,30,15,33,315,630,630,315,33,15,30,1]),
                       dmin=6)
  
  # - K11N30
  SLVD_fer = sim_result(num_sims=np.array([30548, 69853, 183843, 270313, 1e6]),
                   num_errors=np.array([211, 209, 201, 100, 94]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_16_fer = sim_result(num_sims=np.array([30548, 69863, 183843, 582422, 2028557]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_14_fer = sim_result(num_sims=np.array([29990, 69379, 181539, 581798, 2011229]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_12_fer = sim_result(num_sims=np.array([23421, 54558, 136343, 435260, 1503466]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_10_fer = sim_result(num_sims=np.array([15462, 31433, 78293, 214108, 759293]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_8_fer = sim_result(num_sims=np.array([8744, 19097, 42774, 97756, 378923]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_HD_0_fer = sim_result(num_sims=np.array([6693, 14902, 32491, 76692, 272780]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SSD_SLVD_Gabriel_16_fer = sim_result(num_sims=np.array([26922, 61855, 156614, 509878, 1639903]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  N = 30
  K = 11
  # EsNo
  esno_dB = np.arange(-5,2,0.1)
  esno_linear = 10**(esno_dB/10)
  # EbNo
  ebno_linear = esno_linear * N / K
  ebno_dB = 10 * np.log10(ebno_linear)
  # Es/sigma^2
  es_over_sigma_sqrd_linear = esno_linear * 2
  es_over_sigma_sqrd_dB = 10*np.log10(es_over_sigma_sqrd_linear)

  ## - bounds
  dsub_10001 = dsu(ELF_10001, esno_linear)
  dsub_10011 = dsu(CRC_10011, esno_linear)
  dsub_10101 = dsu(ELF_10101, esno_linear)
  dsub_10111 = dsu(ELF_10111, esno_linear)
  dsub_11001 = dsu(CRC_11001, esno_linear)
  dsub_11011 = dsu(ELF_11011, esno_linear)
  dsub_11101 = dsu(ELF_11101, esno_linear)
  dsub_11111 = dsu(CRC_11111, esno_linear)
  
  rcu_filename = 'k11n30_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu = mat_contents['k11n30_rcu']

  # NA
  P_na = []
  for esno in esno_linear:
    p_na = normal_approx(esno, N, K/N)
    P_na.append(p_na)

  
  fig, ax = plt.subplots(figsize=(7, 5.5))

  # - union bounds
  dsu1, = plt.semilogy(ebno_dB, dsub_11001, linewidth=1.5, marker='o', markerfacecolor='none', markevery=5, label=r'$p(x)=31$')
  dsu2, = plt.semilogy(ebno_dB, dsub_11111, linewidth=1.5, marker='s', markerfacecolor='none', markevery=5, label=r'$p(x)=37$')
  dsu3, = plt.semilogy(ebno_dB, dsub_10001, linewidth=1.5, marker='D', markerfacecolor='none', markevery=5, label=r'$p(x)=21$')
  dsu4, = plt.semilogy(ebno_dB, dsub_11101, linewidth=1.5, marker='d', markerfacecolor='none', markevery=5, label=r'$p(x)=35$')
  dsu5, = plt.semilogy(ebno_dB, dsub_10111, linewidth=1.5, marker='p', markerfacecolor='none', markevery=5, label=r'$p(x)=27$')
  dsu6, = plt.semilogy(ebno_dB, dsub_10011, linewidth=1.5, marker='*', markerfacecolor='none', markevery=5, label=r'$p(x)=23$')
  dsu7, = plt.semilogy(ebno_dB, dsub_11011, linewidth=1.5, marker='h', markerfacecolor='none', markevery=5, label=r'$p(x)=33$')
  dsu8, = plt.semilogy(ebno_dB, dsub_10101, linewidth=1.5, marker='8', markerfacecolor='none', markevery=5, label=r'$p(x)=25$')
  dsu_legend = ax.legend(handles=[dsu1, dsu2, dsu3, dsu4, dsu5, dsu6, dsu7, dsu8], loc='lower left', fontsize=12)
  ax.add_artist(dsu_legend)
  
  # - random coding bound
  rcu, = plt.semilogy(rcu[:,1], rcu[:,3], color='k', linestyle='-.', linewidth=1.5, label='Random Coding Union Bound')

  # - normal approximation
  na, = plt.semilogy(ebno_dB, P_na, linestyle='--', linewidth=1.5, label='Normal Approximation')
  ax.legend(handles=[rcu, na], loc='upper right', fontsize=12)

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([2, 6])
  plt.ylim([1e-6, 1e0])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'Probability of codeword error, $P_{cw}$', fontsize=15)
  # plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()

  plt.figure(figsize=(7, 5.5))
  # - sims
  plt.semilogy(SSD_SLVD_HD_0_fer.ebno_dB, SSD_SLVD_HD_0_fer.num_errors/SSD_SLVD_HD_0_fer.num_sims, '-o', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=0$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_HD_8_fer.ebno_dB, SSD_SLVD_HD_8_fer.num_errors/SSD_SLVD_HD_8_fer.num_sims, '-s', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=8$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_HD_10_fer.ebno_dB, SSD_SLVD_HD_10_fer.num_errors/SSD_SLVD_HD_10_fer.num_sims, '-D', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=10$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_HD_12_fer.ebno_dB, SSD_SLVD_HD_12_fer.num_errors/SSD_SLVD_HD_12_fer.num_sims, '-d', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=12$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_HD_14_fer.ebno_dB, SSD_SLVD_HD_14_fer.num_errors/SSD_SLVD_HD_14_fer.num_sims, '-p', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=14$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_HD_16_fer.ebno_dB, SSD_SLVD_HD_16_fer.num_errors/SSD_SLVD_HD_16_fer.num_sims, '-*', 
             linewidth=1, 
             label=r'SSD-SLVD($d_{H}=16$)', 
             markerfacecolor='none')
  plt.semilogy(SLVD_fer.ebno_dB, SLVD_fer.num_errors/SLVD_fer.num_sims, '-', 
             linewidth=1, 
             color='k',
             label=r'SLVD($L=\infty$)', 
             markerfacecolor='none')
  plt.semilogy(SSD_SLVD_Gabriel_16_fer.ebno_dB, SSD_SLVD_Gabriel_16_fer.num_errors/SSD_SLVD_Gabriel_16_fer.num_sims, '-h', 
             linewidth=1, 
             label=r'SSD-SLVD(All Gabriel Neighbors)', 
             markerfacecolor='none')
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([3.5, 5.5])
  plt.ylim([1e-6, 1e0])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel('Frame Error Rate', fontsize=15)
  plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()



  return


if __name__ == "__main__":
  main()