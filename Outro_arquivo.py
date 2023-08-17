import numpy as np
import matplotlib.pyplot as plt

print("Alo")

x = np.arange(0, 10, 0.01)
y = np.exp(np.sin(x))

plt.plot(x, y)
plt.show()