import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
# Metropolis-Hastings algorithm for sampling arbitrary distributions.
# Here we sample sin^2(r)/r^3
# Nice explanation here: http://flynnmichael.com/2015/06/01/my-favorite-algorithm-metropolis-hastings/

samples = np.zeros([10**6+1,2])
currentSample = np.array([0.00001,0.00001])
j = 0
s = 2 # autocorrelation prevention.
while (j < s*10**6+1):
    # Choose random uniform direction.
    direction = np.random.uniform(0,2*np.pi)
    # Choose a distance d with probability exp(-d)
    distance = np.random.exponential(2)
    # Current point
    x1 = currentSample[0]
    y1 = currentSample[1]
    # Compute next point
    x2 = x1 + distance*np.cos(direction)
    y2 = y1 + distance*np.sin(direction)
    # Acceptance criteria: Min(f(x1)/f(x2),1)
    accept = min((np.sin(np.sqrt(x2**2+y2**2))**2*(x1**2+y1**2)**(1.5))/(np.sin(np.sqrt(x1**2+y1**2))**2*(x2**2+y2**2)**(1.5)),1)
    if(accept > np.random.uniform()):
        currentSample = np.array([x2,y2])
        j = j + 1
    else:
        currentSample = np.array([x1,y1])
    # Choose every s steps to prevent autocorrelation.
    if(j%s==0):
        samples[j/s]=currentSample

fig, ax = plt.subplots()
ax.plot(samples[1:,0],samples[1:,1], ',')
ax.set_aspect(1)#'equal', 'datalim')
ax.set_xlim([-6*np.pi,6*np.pi])
ax.set_ylim([-6*np.pi,6*np.pi])