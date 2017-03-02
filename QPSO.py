# This is a standard Quantum Particle Swarm Optimization algorithm

import numpy as np
import matplotlib.pyplot as plt
from keras.models import load_model

decoder = load_model('Compiled/simpleAutoencoders_decoder.h5')
nuclei_data = np.genfromtxt('CRCHistoDataSets/Detection/nucleis_data.dat', delimiter=',')
nuclei_data = nuclei_data.astype('float32') / 255.
nuclei = nuclei_data[4].copy()

plt.imshow(nuclei.reshape(17, 17, 3, order="F"))
plt.title('Original Input Image')
plt.show()

def cost(arg):
    decoded = decoder.predict(arg.reshape(1, 64))
    return np.sqrt(np.sum((decoded - nuclei)**2))

particleNum = 200
maxIteration = 1000
dim = 64
rangeL = 0
rangeR = 15
rangMax = 15

x = (rangeR - rangeL) * np.random.rand(particleNum, dim) + rangeL
x = x.astype('float32')

pbest = x.copy()
gbest = np.zeros((dim), 'float32')
mbest = np.zeros((dim), 'float32')

f_x = np.zeros((particleNum), 'float32')
f_pbest = np.zeros((particleNum), 'float32')

for i in range(0, particleNum, 1):
    f_x[i] = cost(x[i])
    f_pbest[i] = f_x[i]

g = np.argmin(f_pbest)
gbest = pbest[g]
f_gbest = f_pbest[g]

MINIMUM = f_gbest

for t in range(0, maxIteration, 1):

    beta = 0.5 * (maxIteration - t) / maxIteration + 0.5
    mbest = np.sum(pbest, axis = 0) / particleNum

    for i in range(0, particleNum, 1):

        fi = np.random.rand(dim)
        p = fi * pbest[i] + (1 - fi) * gbest
        u = np.random.rand(dim)
        b = beta * np.absolute(mbest - x[i])
        v = -1 * np.log(u)
        y = p + ((-1) ** np.ceil(0.5 + np.random.rand(dim))) * b * v
        x[i] = np.sign(y) * np.minimum(np.absolute(y), rangMax)
        x[i] = np.absolute(x[i])
        f_x[i] = cost(x[i])

        if f_x[i] < f_pbest[i]:
            pbest[i] = x[i]
            f_pbest[i] = f_x[i]
        if f_pbest[i] < f_gbest:
            gbest = pbest[i]
            f_gbest = f_pbest[i]

        MINIMUM = f_gbest

    print(str(t) + " " + str(MINIMUM))
    to_show = decoder.predict(gbest.reshape(1, 64))
    plt.imshow(to_show.reshape(17, 17, 3, order="F"))
    plt.title('Loss: ' + str(MINIMUM))
    plt.savefig('Animation' + '\img' + str(t) + '.png')
    plt.clf()
    plt.cla()
    plt.close('all')

print(MINIMUM)
print(gbest)
to_show = decoder.predict(gbest.reshape(1, 64))
plt.imshow(to_show.reshape(17, 17, 3, order="F"))
plt.show()







