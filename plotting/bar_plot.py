import matplotlib.pyplot as plt
import numpy as np
flops = np.array([19360, 193450])
planeTouch = np.array([13149, 7367, 13581])
planeExplode = np.array([412810, 130015, 188922])

perfTouch = flops[0] * 1/planeTouch
perfExplode = flops[1] * 1/planeExplode

data_perf = np.array([perfTouch, perfExplode])
data_cycles = np.array([planeTouch, planeExplode])

# langs = ['Baseline', 'Vec', 'F16C']
# plt.bar(langs, data_perf[0], color = 'g', width = 0.25)
# plt.bar(langs, data_perf[1], color = 'r', width = 0.25)

ind = np.arange(3)  # the x locations for the groups
width = 0.5

fig = plt.figure()
ax = fig.add_subplot(111)
langs = ['Baseline', 'Vec', 'F16C']
ax.bar(ind, data_perf[0], color = 'g', width = 0.8*width)
ax.bar(ind + width, data_perf[1], color = 'r', width = 0.8*width)
ax.set_xticks(ind + width / 2)
ax.set_xticklabels( ('Baseline', 'Vec', 'F16C') )


font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 9
}

plt.xlabel('Operational Intensity [Flops / Byte]', fontdict=font)
plt.grid(axis='y')
plt.title('Intel® Xeon D-2141i @ 2.2 GHz'
          '\nL1: 512KiB L2:8MiB L3:11MiB'
          '\nMemory Bandwidth: 30.9 bytes/cycle'
          '\nCompiler: GCC 11.1'
          '\nCompiler Flags: -O3 -march=skylake-avx512 -ffast-math'
          '\n\n'
          'Performance [flops/cycle]',
          # 'Cycles',
          fontdict=font, loc='left')

plt.legend(labels=['PlaneBodyTouchingWings', 'PlaneBodyWingsExploded'], bbox_to_anchor=(1.05, 1))
plt.show()
