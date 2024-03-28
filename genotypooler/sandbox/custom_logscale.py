import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import FormatStrFormatter
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms

class CustomScale(mscale.ScaleBase):
    name = 'custom'

    def __init__(self, axis='x'):  #, axis, **kwargs):
        mscale.ScaleBase.__init__(self, axis)
        self.thresh = None #thresh

    def get_transform(self):
        return self.CustomTransform(self.thresh)

    def set_default_locators_and_formatters(self, axis):
        pass

    class CustomTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            # return np.log(1+a)
            return np.log10(1 - a)

        def inverted(self):
            return CustomScale.InvertedCustomTransform(self.thresh)

    class InvertedCustomTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            # return np.exp(a)-1
            return 1 - np.power(10, a)

        def inverted(self):
            return CustomScale.CustomTransform(self.thresh)

# Now that the Scale class has been defined, it must be registered so
# that ``matplotlib`` can find it.
mscale.register_scale(CustomScale)

u = np.arange(0.05, 0.55, 0.5 / 5)
z = np.arange(0.1, 5.1, 1.0)
print(u)
thick = 1 - np.power(10, -1 * z[::-1])
thick = [0.999, 0.99, 0.9, 0.8, 0.4]
#thick = [0.001, 0.01, 0.1, 0.2, 0.6]
print(thick)

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(111)
ax1.plot(u, thick, marker='o', linewidth=2, c='k')

plt.xlabel(r'$\rm{redshift}$', size=16)
plt.ylabel(r'$\rm{thickness\ (kpc)}$', size=16)
plt.gca().set_yscale('custom')
plt.show()