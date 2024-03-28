
# Implementation of matplotlib function
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(range(10))
axin1 = ax.inset_axes([0.8, 0.1,
                       0.15, 0.15])
axin1.plot(range(2), c='g')

axin2 = ax.inset_axes([5, 7, 2.3, 2.3], transform=ax.transData)
axin2.plot(range(2), c='tab:orange')

ax.set_title('matplotlib.axes.Axes.inset_axes() Example',
             fontsize = 14, fontweight ='bold')
plt.show()
