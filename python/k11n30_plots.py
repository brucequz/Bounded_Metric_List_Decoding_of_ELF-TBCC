from dataclasses import dataclass
from bounds import dsu, normal_approx
from matplotlib.patches import Ellipse
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# 1. Fix the Font Type (Type 42 is TrueType, avoids the Type 3 error)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# 2. Set font family to Sans-Serif (like MATLAB's Helvetica/Arial)
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans']


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
  max_list_size: Optional[np.array] = None
  avg_list_size: Optional[np.array] = None




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
  # SLVD_fer = sim_result(num_sims=np.array([30548, 69853, 183843, 270313, 1e6]),
  #                  num_errors=np.array([211, 209, 201, 100, 94]),
  #                  ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  SLVD_fer = sim_result(num_sims=np.array([15091, 30487, 76864, 175821, 642491, 2093865, 7870474]),
                   num_errors=np.array([200, 200, 200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3, 3.5, 4, 4.5, 5, 5.5, 6]),
                   max_list_size=np.array([3445, 3159, 2726, 2948, 2150, 1550, 1004]),
                   avg_list_size=np.array([16.61911, 9.062486, 4.662989, 2.5275877, 1.6250532, 1.2495171, 1.100459]))
  
  TB_SSD_SLVD_HD_16_fer = sim_result(num_sims=np.array([30548, 69863, 183843, 582422, 2028557]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_16_fer = sim_result(num_sims=np.array([30487, 76864, 183843, 653068, 2039663]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_HD_14_fer = sim_result(num_sims=np.array([29990, 69379, 181539, 581798, 2011229]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_14_fer = sim_result(num_sims=np.array([30403, 75700, 172940, 630573, 2087067]),
                  num_errors=np.array([200, 200, 200, 200, 200]),
                  ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_HD_12_fer = sim_result(num_sims=np.array([23421, 54558, 136343, 435260, 1503466]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_12_fer = sim_result(num_sims=np.array([22042, 51088, 107732, 344975, 1048015]),
                  num_errors=np.array([200, 200, 200, 200, 200]),
                  ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_HD_10_fer = sim_result(num_sims=np.array([15462, 31433, 78293, 214108, 759293]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_10_fer = sim_result(num_sims=np.array([9888, 17572, 33386, 68097, 187825]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_HD_8_fer = sim_result(num_sims=np.array([8744, 16859, 42224, 104498, 412492]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_8_fer = sim_result(num_sims=np.array([4797, 7617, 14769, 26679, 73818]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_HD_0_fer = sim_result(num_sims=np.array([6693, 13245, 34984, 81429, 286258]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  nonTB_SSD_SLVD_HD_0_fer = sim_result(num_sims=np.array([2964, 5050, 8692, 17433, 38991]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  TB_SSD_SLVD_Gabriel_16_fer = sim_result(num_sims=np.array([26922, 61855, 156614, 509878, 1639903]),
                   num_errors=np.array([200, 200, 200, 200, 200]),
                   ebno_dB=np.array([3.5, 4, 4.5, 5, 5.5]))
  
  N = 30
  K = 11
  R = K/N
  # EbNo
  ebno_dB = np.arange(1, 7.1, 0.1)
  ebno_linear = 10**(0.1*ebno_dB)
  # EsNo
  esno_linear = ebno_linear * R
  esno_dB = 10*np.log10(esno_linear)
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

  random_coding_k11n30 = sim_result(num_sims=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
              num_errors=np.array([0.1285, 0.0834, 5.28e-2, 2.87e-2, 1.61e-2, 8.34e-3, 3.37e-3, 1.47e-3, 6.19e-4, 2.45e-04, 9.01e-5]),
              ebno_dB=np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]))
  


  # NA
  P_na = []
  for esno in esno_linear:
    p_na = normal_approx(esno, N, K/N)
    P_na.append(p_na)

  # np.savez("k11n30_10011_dsub.npz", ebno_dB=ebno_dB, esno_linear=esno_linear, K=K, N=N, dsub=dsub_10011)
  
  fig, ax = plt.subplots(figsize=(7, 5))
  ax.set_yscale('log')
  # - union bounds
  dsu1, = plt.semilogy(ebno_dB, dsub_11001, linewidth=1.5, marker='o', markerfacecolor='none', markevery=5, label=r'$g_e(x)=31, A_{6}=45$')
  dsu2, = plt.semilogy(ebno_dB, dsub_11111, linewidth=1.5, marker='s', markerfacecolor='none', markevery=5, label=r'$g_e(x)=37, A_{6}=30$')
  dsu3, = plt.semilogy(ebno_dB, dsub_10001, linewidth=1.5, marker='D', markerfacecolor='none', markevery=5, label=r'$g_e(x)=21, A_{6}=4$')
  dsu4, = plt.semilogy(ebno_dB, dsub_11101, linewidth=1.5, marker='d', markerfacecolor='none', markevery=5, label=r'$g_e(x)=35, A_{6}=4$')
  dsu5, = plt.semilogy(ebno_dB, dsub_10111, linewidth=1.5, marker='p', markerfacecolor='none', markevery=5, label=r'$g_e(x)=27, A_{6}=3$')
  dsu6, = plt.semilogy(ebno_dB, dsub_10011, linewidth=1.5, marker='*', markerfacecolor='none', markevery=5, label=r'$g_e(x)=23, A_{8}=30$')
  dsu7, = plt.semilogy(ebno_dB, dsub_11011, linewidth=1.5, marker='v', markerfacecolor='none', markevery=5, label=r'$g_e(x)=33, A_{8}=21$')
  dsu8, = plt.semilogy(ebno_dB, dsub_10101, linewidth=1.5, marker='^', markerfacecolor='none', markevery=5, label=r'$g_e(x)=25, A_{8}=20$')
  dsu_legend = ax.legend(handles=[dsu1, dsu2, dsu3, dsu4, dsu5, dsu6, dsu7, dsu8], loc='lower left', fontsize=12)
  ax.add_artist(dsu_legend)
  
  # - random coding bound
  # rcu, = plt.semilogy(rcu[:,1], rcu[:,3], color='k', linestyle='-.', linewidth=1.5, label='Random Coding Union Bound')

  # rcu_sim, = plt.semilogy(random_coding_k11n30.ebno_dB, random_coding_k11n30.num_errors/random_coding_k11n30.num_sims, ':s', 
        # linewidth=1,
        # color='k',
        # label=r'Random Coding Simulation')

  # - normal approximation
  # na, = plt.semilogy(ebno_dB, P_na, linestyle='--', linewidth=1.5, label='Normal Approximation')
  # ax.legend(handles=[rcu_sim, rcu, na], loc='upper right', fontsize=12)

  # ellipses and arrows
  plt.annotate(
    r'$31$',      # The text
    xy=(6.52, 2.1e-4),       # Tip of the arrow (x, y)
    xytext=(6.68, 2.1e-4),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='left'           # Align text to the right of the xytext point
  )

  plt.annotate(
    r'$37$',      # The text
    xy=(6.48, 1.4e-4),       # Tip of the arrow (x, y)
    xytext=(6.35, 1.4e-4),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='right'           # Align text to the right of the xytext point
  )

  cx, cy = 6.5, 2e-5
  h = cy * (10**0.10 - 10**-0.10)
  ax.add_patch(Ellipse((cx, cy), 0.4, h, edgecolor='k', fc='none', lw=2))
  ax.text(6.8, 2.5e-5, "21, 35, 27", color='black', fontsize=11, ha='center')

  plt.annotate(
    r'$23$',      # The text
    xy=(6.52, 5e-6),       # Tip of the arrow (x, y)
    xytext=(6.68, 5e-6),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='left'           # Align text to the right of the xytext point
  )

  plt.annotate(
    r'$33,25$',      # The text
    xy=(6.48, 3.8e-6),       # Tip of the arrow (x, y)
    xytext=(6.35, 3.8e-6),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='right'           # Align text to the right of the xytext point
  )

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([4, 7])
  plt.ylim([1e-6, 1e-2])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'Probability of codeword error, $P_{cw}$', fontsize=15)
  # plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()

  fig, ax = plt.subplots(figsize=(7, 4.5))
  # - sims
  ssd_hd0, = plt.semilogy(TB_SSD_SLVD_HD_0_fer.ebno_dB, TB_SSD_SLVD_HD_0_fer.num_errors/TB_SSD_SLVD_HD_0_fer.num_sims, '-o', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 0)$', 
             markerfacecolor='none')
  ssd_hd8, = plt.semilogy(TB_SSD_SLVD_HD_8_fer.ebno_dB, TB_SSD_SLVD_HD_8_fer.num_errors/TB_SSD_SLVD_HD_8_fer.num_sims, '-s', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 8)$', 
             markerfacecolor='none')
  ssd_hd10, = plt.semilogy(TB_SSD_SLVD_HD_10_fer.ebno_dB, TB_SSD_SLVD_HD_10_fer.num_errors/TB_SSD_SLVD_HD_10_fer.num_sims, '-D', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 10)$', 
             markerfacecolor='none')
  ssd_hd12, = plt.semilogy(TB_SSD_SLVD_HD_12_fer.ebno_dB, TB_SSD_SLVD_HD_12_fer.num_errors/TB_SSD_SLVD_HD_12_fer.num_sims, '-d', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 12)$', 
             markerfacecolor='none')
  ssd_hd14, = plt.semilogy(TB_SSD_SLVD_HD_14_fer.ebno_dB, TB_SSD_SLVD_HD_14_fer.num_errors/TB_SSD_SLVD_HD_14_fer.num_sims, '-p', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 14)$', 
             markerfacecolor='none')
  ssd_hd16, = plt.semilogy(TB_SSD_SLVD_HD_16_fer.ebno_dB, TB_SSD_SLVD_HD_16_fer.num_errors/TB_SSD_SLVD_HD_16_fer.num_sims, '-*', 
             linewidth=1, 
             label=r'SSD ${\cal N} (W_\text{H} \leq 16)$', 
             markerfacecolor='none')
  slvd, = plt.semilogy(SLVD_fer.ebno_dB, SLVD_fer.num_errors/SLVD_fer.num_sims, '-', 
             linewidth=2, 
             color='k',
             label=r'SLVD (ML Performance)', 
             markerfacecolor='none')
  ssd_gn, = plt.semilogy(TB_SSD_SLVD_Gabriel_16_fer.ebno_dB, TB_SSD_SLVD_Gabriel_16_fer.num_errors/TB_SSD_SLVD_Gabriel_16_fer.num_sims, '--h', 
             linewidth=1, 
             label=r'SSD ${\cal N} = {\cal G}$', 
             markerfacecolor='none')
  ssd_legend = ax.legend(ncol=2, handles=[ssd_hd0, ssd_hd8, ssd_hd10, ssd_hd12, ssd_hd14, ssd_hd16], loc='upper right', fontsize=12)
  ax.add_artist(ssd_legend)
  ax.legend(handles=[ssd_gn, slvd], loc='lower left', fontsize=12)

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([3.5, 5.5])
  plt.ylim([5e-5, 1e-1])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'Probability of codeword error, $P_{cw}$', fontsize=15)
  plt.tight_layout()
  plt.show()

  fig, ax = plt.subplots(figsize=(7, 4.5))
  # - sims

  nontb_ssd_hd0, = plt.semilogy(nonTB_SSD_SLVD_HD_0_fer.ebno_dB, nonTB_SSD_SLVD_HD_0_fer.num_errors/nonTB_SSD_SLVD_HD_0_fer.num_sims, '-o', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 0)$', 
             markerfacecolor='none')
  
  nontb_ssd_hd8, = plt.semilogy(nonTB_SSD_SLVD_HD_8_fer.ebno_dB, nonTB_SSD_SLVD_HD_8_fer.num_errors/nonTB_SSD_SLVD_HD_8_fer.num_sims, '-s', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 8)$', 
             markerfacecolor='none')
  
  nontb_ssd_hd10, = plt.semilogy(nonTB_SSD_SLVD_HD_10_fer.ebno_dB, nonTB_SSD_SLVD_HD_10_fer.num_errors/nonTB_SSD_SLVD_HD_10_fer.num_sims, '-D', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 10)$', 
             markerfacecolor='none')
  
  nontb_ssd_hd12, = plt.semilogy(nonTB_SSD_SLVD_HD_12_fer.ebno_dB, nonTB_SSD_SLVD_HD_12_fer.num_errors/nonTB_SSD_SLVD_HD_12_fer.num_sims, '-d', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 12)$', 
             markerfacecolor='none')
  
  nontb_ssd_hd14, = plt.semilogy(nonTB_SSD_SLVD_HD_14_fer.ebno_dB, nonTB_SSD_SLVD_HD_14_fer.num_errors/nonTB_SSD_SLVD_HD_14_fer.num_sims, '-p', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 14)$', 
             markerfacecolor='none')
  
  nontb_ssd_hd16, = plt.semilogy(nonTB_SSD_SLVD_HD_16_fer.ebno_dB, nonTB_SSD_SLVD_HD_16_fer.num_errors/nonTB_SSD_SLVD_HD_16_fer.num_sims, '-*', 
             linewidth=1, 
             label=r'nonTB SSD ${\cal N} (W_\text{H} \leq 16)$', 
             markerfacecolor='none')
  
  slvd, = plt.semilogy(SLVD_fer.ebno_dB, SLVD_fer.num_errors/SLVD_fer.num_sims, '-', 
            linewidth=2, 
            color='k',
            label=r'SLVD (ML Performance)', 
            markerfacecolor='none')

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([3.5, 5.5])
  plt.ylim([5e-5, 1e-1])
  plt.legend(fontsize=12, ncol=2, loc='lower left')
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'Probability of codeword error, $P_{cw}$', fontsize=15)
  plt.tight_layout()
  plt.show()


  fig, ax = plt.subplots(figsize=(7, 4.5))
  ax.set_yscale('log')
  plt.semilogy(SLVD_fer.ebno_dB, SLVD_fer.max_list_size, '.', linestyle='-', lw=1.5)
  plt.semilogy(SLVD_fer.ebno_dB, SLVD_fer.avg_list_size, '.', linestyle='-', lw=1.5)

  plt.annotate(
    r'$\text{Max List Size}$',      # The text
    xy=(4.65, 3e3),       # Tip of the arrow (x, y)
    xytext=(5.0, 4e3),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='bottom',         # Center text vertically with arrow
    ha='center',           # Align text to the right of the xytext point
  )

  plt.annotate(
    r'$\text{Average List Size}$',      # The text
    xy=(4.8, 2e0),       # Tip of the arrow (x, y)
    xytext=(5.0, 4e0),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1.5           # Line width
    ),
    va='bottom',         # Center text vertically with arrow
    ha='center',           # Align text to the right of the xytext point
  )

  plt.xlim([3, 6])
  plt.ylim([1e-0, 1e+4])
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'List Size', fontsize=15)
  plt.tight_layout()
  plt.show()

  return


if __name__ == "__main__":
  main()