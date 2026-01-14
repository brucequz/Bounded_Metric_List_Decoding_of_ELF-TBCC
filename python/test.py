import numpy as np
import matplotlib.pyplot as plt
from bounds import bi_awgn_capacity_stable, normal_approx

def main():
  
  # ebno_dB = np.arange(-3, 8, 0.1)
  # ebno_linear = 10**(0.1*ebno_dB)
  # print(ebno_linear)

  

  ebno_linear = []
  esno_dB = np.arange(-5, 8, 0.1)
  esno_linear = 10**(0.1*esno_dB)
  rho_linear  = 2 * esno_linear
  capacity = []

  for r in rho_linear:
    print(r)
    C = bi_awgn_capacity_stable(r)
    R = C
    print("C: ", C)
    capacity.append(C)
    ebno_linear.append(1/(2*R)*r)
  
  ebno_dB = 10*np.log10(np.array(ebno_linear))
  
  plt.figure()
  plt.plot(ebno_dB, capacity, label='soft')
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([-2,10])
  plt.ylim([0, 1])
  plt.yticks(np.arange(0,1.1,0.1))
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel('Capacity', fontsize=15)
  plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()

  k=12
  n=24
  R=0.5
  ebno_dB=np.arange(1,5.1,0.1)
  ebno_linear=10**(0.1*ebno_dB)

  esno_linear = ebno_linear * R
  P_na = []
  for esno in esno_linear:
    p_na = normal_approx(esno, n, R)
    P_na.append(p_na)
    
  plt.figure()
  plt.semilogy(ebno_dB, P_na, label='Normal approximation')
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)
  plt.xlim([1, 5])
  plt.ylim([1e-6, 1e-0])
  plt.xlabel(r'$\frac{E_b}{N_o} (\mathrm{dB})$', fontsize=15)
  plt.ylabel(r'$P_e$', fontsize=15)
  plt.legend(fontsize=12)
  plt.tight_layout()
  plt.show()





  return

if __name__ == "__main__":
  main()