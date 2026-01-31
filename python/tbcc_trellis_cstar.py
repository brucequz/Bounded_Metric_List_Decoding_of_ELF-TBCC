import numpy as np
import matplotlib.pyplot as plt


def main():
  TrellisCstar_to_received_vec = np.fromfile("../output/K15N30/trellis_cstar_to_received.txt", sep=" ")
  TBCstar_to_received_vec = np.fromfile("../output/K15N30/tbcc_cstar_to_received.txt", sep=" ")

  trellis_avg = np.average(TrellisCstar_to_received_vec)
  trellis_var = np.var(TrellisCstar_to_received_vec)
  print("trellis_avg: ", trellis_avg)
  print("trellis_var: ", trellis_var)
  tbcc_avg = np.average(TBCstar_to_received_vec)
  tbcc_var = np.var(TBCstar_to_received_vec)
  print("tbcc_avg: ", tbcc_avg)
  print("tbcc_var: ", tbcc_var)

  plt.figure(figsize=(7, 5.5))
  plt.hist(TrellisCstar_to_received_vec, bins=100, density=True, color='skyblue', edgecolor='black', label="Trellis-based SSD")
  plt.hist(TBCstar_to_received_vec, bins=100, density=True, color='lightgreen', edgecolor='black', label="TBCC-based SSD")
  plt.ylim([0, 1])
  plt.title(r"Distribution of $|c_t^* - y|$")
  plt.xlabel("Euclidean Distance")
  plt.ylabel("Probability")
  plt.grid(axis='y', alpha=0.75)
  plt.legend()
  plt.show()


if __name__ == "__main__":
  main()