import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

a = loadmat('Allcos.mat')

st = [[element for element in upperElement] for upperElement in a['Allcos']]

print(type(st))
print(np.shape(np.array(st)))

st = np.array(st)

st1 = st[0]
st2 = st[1]
st3 = st[2]
st4 = st[3]
st5 = st[4]
plt.title("Cosine angles between queries and documents")
plt.xlabel("Document")
plt.ylabel("Value of cosine")

plt.plot(st1, label='Query 1')
plt.plot(st2, label='Query 2')
plt.plot(st3, label='Query 3')
plt.plot(st4, label='Query 4')
plt.plot(st5, label='Query 5')
plt.legend()

plt.savefig('plot.png')
plt.show()