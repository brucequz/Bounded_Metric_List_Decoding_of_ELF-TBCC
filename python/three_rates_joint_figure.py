import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import matplotlib

# 1. Fix the Font Type (Type 42 is TrueType, avoids the Type 3 error)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# 2. Set font family to Sans-Serif (like MATLAB's Helvetica/Arial)
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans']
import scipy
from dataclasses import dataclass

earthy_palette = {
    'Forest': ['#0A260D', '#1B5E20', '#2E7D32', '#4CAF50', '#81C784', '#C8E6C9'],
    'Slate':  ['#10171B', '#263238', '#455A64', '#607D8B', '#90A4AE', '#CFD8DC'],
    'Terra':  ['#3E2723', '#5D4037', '#8D6E63', '#A1887F', '#D7CCC8', '#EFEBE9']
}

professional_palette = {
    'Indigo': ['#1A237E', '#303F9F', '#3F51B5', '#7986CB', '#9FA8DA', '#C5CAE9'],
    'Teal':   ['#004D40', '#00796B', '#009688', '#4DB6AC', '#80CBC4', '#B2DFDB'],
    'Amber':  ['#E65100', '#F57C00', '#FF9800', '#FFB74D', '#FFCC80', '#FFE0B2']
}

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
  
  # - k11n30
  rcu_filename = 'k11n30_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu_k11n30 = mat_contents['k11n30_rcu']
  na_k11n30 = np.load("k11n30_na.npz")
  SLVD_K11N30_fer = sim_result(num_sims=np.array([15091, 30487, 76864, 175821, 642491, 2093865, 7870474]),
                num_errors=np.array([200, 200, 200, 200, 200, 200, 200]),
                ebno_dB=np.array([3, 3.5, 4, 4.5, 5, 5.5, 6]))
  dsub_10011_k11n30 = np.load("k11n30_10011_dsub.npz")

  

  # - k21n62
  rcu_filename = 'k21n62_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu_k21n62 = mat_contents['k21n62_rcu']
  na_k21n62 = np.load("k21n62_na.npz")
  SLVD_K21N62_fer = sim_result(num_sims=np.array([134732, 506099, 3026802, 18684105, 165563627]),
                              num_errors=np.array([200, 200, 200, 200, 200]),
                              ebno_dB=np.array([3, 3.5, 4, 4.5, 5]))
  dsub_11101101001 = np.load("k21n62_11101101001_dsub.npz")

  # - k51n126
  rcu_filename = 'k51n126_rcu.mat'
  mat_contents = scipy.io.loadmat(rcu_filename)
  rcu_k51n126 = mat_contents['k51n126_rcu']
  na_k51n126 = np.load("k51n126_na.npz")
  SLVD_K51N126_fer = sim_result(num_sims=np.array([52934, 137818, 319086, 991788, 2727906, 9225120, 37343200, 101582207]),
                              num_errors=np.array([200, 200, 200, 200, 200, 200, 200, 200]),
                              ebno_dB=np.array([2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75]))
  
  SLVD_K51N126v8_fer = sim_result(num_sims=np.array([12804, 63200, 159552, 432676, 18693847]),
                                num_errors=np.array([200, 200, 200, 200, 200]),
                                ebno_dB=np.array([1.5, 2, 2.25, 2.5, 3.25]))

  fig, ax = plt.subplots(figsize=(7,4.5))
  ax.set_yscale('log')
  plt.semilogy(rcu_k11n30[:,1], rcu_k11n30[:,3], color='r', linestyle='-.', linewidth=1.5, label='(30, 11) Random Coding Union Bound')
  plt.semilogy(na_k11n30['ebno_dB'], na_k11n30['P_na'], color='b', linestyle='--', linewidth=1.5, label='(30, 11) Normal Approximation')
  plt.semilogy(SLVD_K11N30_fer.ebno_dB, SLVD_K11N30_fer.num_errors/SLVD_K11N30_fer.num_sims, '-^', 
            linewidth=2, 
            color='k',
            label=r'(30, 11) SLVD (ML Performance)', 
            markerfacecolor='none')
  plt.semilogy(dsub_10011_k11n30['ebno_dB'], dsub_10011_k11n30['dsub'], '-', 
            linewidth=1.5, 
            marker='v',
            markevery=5,
            color='k',
            label=r'(30, 11) Distance Spectrum Union Bound', 
            markerfacecolor='none')
  cx, cy = 3.9, 3e-3
  h = cy * (10**0.1 - 10**-0.1)
  ax.add_patch(Ellipse((cx, cy), 0.8, h, edgecolor='k', fc='none', lw=2))
  ax.text(4.7, 2e-3, "(30, 11)", color='black', fontsize=11, ha='center')

  # - k21n62
  plt.semilogy(rcu_k21n62[:,1], rcu_k21n62[:,3], color='r', linestyle='-.', linewidth=1.5, label='(62, 21) Random Coding Union Bound')
  plt.semilogy(na_k21n62['ebno_dB'], na_k21n62['P_na'], color='b', linestyle='--', linewidth=1.5, label='(62, 21) Normal Approximation')
  plt.semilogy(SLVD_K21N62_fer.ebno_dB, SLVD_K21N62_fer.num_errors/SLVD_K21N62_fer.num_sims, '-^', 
            linewidth=2, 
            color='k',
            label=r'(62, 21), SLVD (ML Performance)', 
            markerfacecolor='none')
  plt.semilogy(dsub_11101101001['ebno_dB'], dsub_11101101001['dsub'], '-', 
          linewidth=1.5, 
          color='k',
          marker='v',
          markevery=5,
          label=r'(62, 21) Distance Spectrum Union Bound', 
          markerfacecolor='none')
  cx, cy = 3.85, 1.5e-4
  h = cy * (10**0.1 - 10**-0.1)
  ax.add_patch(Ellipse((cx, cy), 0.4, h, edgecolor='k', fc='none', lw=2))
  ax.text(3.5, 7e-5, "(62, 21)", color='black', fontsize=11, ha='center')

  # -k51n126
  plt.semilogy(rcu_k51n126[:,1], rcu_k51n126[:,3], color='r', linestyle='-.', linewidth=1.5, label='(126, 51) Random Coding Union Bound')
  plt.semilogy(na_k51n126['ebno_dB'], na_k51n126['P_na'], color='b', linestyle='--', linewidth=1.5, label='(126, 51) Normal Approximation')
  plt.semilogy(SLVD_K51N126_fer.ebno_dB, SLVD_K51N126_fer.num_errors/SLVD_K51N126_fer.num_sims, '-^', 
          linewidth=2, 
          color='k',
          label=r'(126, 51), $\nu=6$, SLVD (ML Performance)', 
          markerfacecolor='none')
  plt.annotate(
    r'$\nu = 6$',      # The text
    xy=(3.55, 4e-6),       # Tip of the arrow (x, y)
    xytext=(3.8, 4e-6),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='left'           # Align text to the right of the xytext point
  )
  plt.semilogy(SLVD_K51N126v8_fer.ebno_dB, SLVD_K51N126v8_fer.num_errors/SLVD_K51N126v8_fer.num_sims, '-^', 
          linewidth=2, 
          color='k',
          label=r'(126, 51), $\nu=8$, SLVD (ML Performance)', 
          markerfacecolor='none')
  plt.annotate(
    r'$\nu = 8$',      # The text
    xy=(3, 3e-5),       # Tip of the arrow (x, y)
    xytext=(2.75, 3e-5),   # Start of the text (x, y)
    arrowprops=dict(
        arrowstyle='->', # Arrow shape
        color='black',     # Arrow color
        lw=1           # Line width
    ),
    va='center',         # Center text vertically with arrow
    ha='right'           # Align text to the right of the xytext point
  )
  cx, cy = 2.35, 1e-3
  h = cy * (10**0.1 - 10**-0.1)
  ax.add_patch(Ellipse((cx, cy), 0.4, h, edgecolor='k', fc='none', lw=2))
  ax.text(2, 4.5e-4, "(126, 51)", color='black', fontsize=11, ha='center')

  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([1, 7])
  plt.ylim([1e-6, 1e-1])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'Probability of codeword error, $P_{cw}$', fontsize=15)
  legend_elements = [
    Line2D([0], [0], color='r', lw=2, linestyle='-.', label='Random Coding Union Bound'),
    Line2D([0], [0], color='black', lw=2, linestyle='-', marker='v', markerfacecolor='none', label='Distance Spectrum Union Bound'),
    Line2D([0], [0], color='black', lw=2, linestyle='-', marker='^', markerfacecolor='none', label='Simulation (ML Performance)'),
    Line2D([0], [0], color='b', lw=2, linestyle='--', label='Normal Approximation')
  ]
  plt.legend(handles=legend_elements, loc='upper right')
  plt.tight_layout()
  plt.show()



  return

if __name__ == "__main__":
  main()