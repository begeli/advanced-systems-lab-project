import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



def calculateSpeedup(optimized, baseline):
    speedup = []
    for i in range(0, len(baseline)):
        speedup.append(baseline[i] / optimized[i])

    return speedup

def calculatePerf(cycles, flops):
    perf = []
    for i in range(4, 23, 1):
        perf.append(flops[i - 4] / cycles[i - 4])

    return perf

def calculateOI(flops):
    op_intensity = []

    for j in range(4,23,1):
        bytesTransferred = 2 * 12 * 2**j if j <= 18 else 4 * 2 * 12 * 2**j
        intensity = flops[j - 4] / bytesTransferred
        op_intensity.append(intensity)

    return op_intensity





basehover = pd.read_csv('../Data/bvh/BASEtoyotahover-multires-data.csv', delimiter=';')
vechover = pd.read_csv('../Data/bvh/VECtoyotahover-multires-data.csv', delimiter=';')
f16chover = pd.read_csv('../Data/bvh/F16Ctoyotahover-multires-data.csv', delimiter=';')
basehover = pd.read_csv('../Data/bvh/BASEtoyotahover-multires-data.csv', delimiter=';')
vechover = pd.read_csv('../Data/bvh/VECtoyotahover-multires-data.csv', delimiter=';')
f16chover = pd.read_csv('../Data/bvh/F16Ctoyotahover-multires-data.csv', delimiter=';')

basetouch = pd.read_csv('../Data/bvh/BASEtoyotatouch-multires-data.csv', delimiter=';')
vectouch = pd.read_csv('../Data/bvh/VECtoyotatouch-multires-data.csv', delimiter=';')
f16ctouch = pd.read_csv('../Data/bvh/F16Ctoyotatouch-multires-data.csv', delimiter=';')
basetouch = pd.read_csv('../Data/bvh/BASEtoyotatouch-multires-data.csv', delimiter=';')
vectouch = pd.read_csv('../Data/bvh/VECtoyotatouch-multires-data.csv', delimiter=';')
f16ctouch = pd.read_csv('../Data/bvh/F16Ctoyotatouch-multires-data.csv', delimiter=';')

print(basehover)

triangles = np.flip(basehover['triangles'].to_numpy())
print(triangles)

basehovercycles = np.flip(basehover['cycles'].to_numpy())
basehoverflops = np.flip(basehover['flops'].to_numpy())
basehovermemory = np.flip(basehover['memory'].to_numpy())
basehovertriangletests = np.flip(basehover['triangletests'].to_numpy())

vechovercycles = np.flip(vechover['cycles'].to_numpy())
vechoverflops = np.flip(vechover['flops'].to_numpy())
vechovermemory = np.flip(vechover['memory'].to_numpy())
vechovertriangletests = np.flip(vechover['triangletests'].to_numpy())

f16chovercycles = np.flip(f16chover['cycles'].to_numpy())
f16choverflops = np.flip(f16chover['flops'].to_numpy())
f16chovermemory = np.flip(f16chover['memory'].to_numpy())
f16chovertriangletests = np.flip(f16chover['triangletests'].to_numpy())


vechoverspeedup = calculateSpeedup(vechovercycles, basehovercycles)
f16choverspeedup = calculateSpeedup(f16chovercycles, basehovercycles)



basetouchcycles = np.flip(basetouch['cycles'].to_numpy())
basetouchflops = np.flip(basetouch['flops'].to_numpy())
basetouchmemory = np.flip(basetouch['memory'].to_numpy())
basetouchtriangletests = np.flip(basetouch['triangletests'].to_numpy())

vectouchcycles = np.flip(vectouch['cycles'].to_numpy())
vectouchflops = np.flip(vectouch['flops'].to_numpy())
vectouchmemory = np.flip(vectouch['memory'].to_numpy())
vectouchtriangletests = np.flip(vectouch['triangletests'].to_numpy())

f16ctouchcycles = np.flip(f16ctouch['cycles'].to_numpy())
f16ctouchflops = np.flip(f16ctouch['flops'].to_numpy())
f16ctouchmemory = np.flip(f16ctouch['memory'].to_numpy())
f16ctouchtriangletests = np.flip(f16ctouch['triangletests'].to_numpy())


vectouchspeedup = calculateSpeedup(vectouchcycles, basetouchcycles)
f16ctouchspeedup = calculateSpeedup(f16ctouchcycles, basetouchcycles)



plt.plot(triangles, basehovercycles, 'o--', label='Baseline (disjoint)', color='tab:blue')
plt.plot(triangles, basetouchcycles, 'o-', label='Baseline (intersecting)', color='tab:blue')
plt.plot(triangles, vechovercycles, 'x--', label='Vectorized (disjoint)', markersize=10, color='tab:orange')
plt.plot(triangles, vectouchcycles, 'x-', label='Vectorized (intersecting)', markersize=10, color='tab:orange')
plt.plot(triangles, f16chovercycles, '^--', label='Half-Precision (disjoint)', color='tab:green')
plt.plot(triangles, f16ctouchcycles, '^-', label='Half-Precision (intersecting)', color='tab:green')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
a_string = "BVH Runtime in Cycles vs. Input Size "
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Runtime in Cycles vs. Input Size "
bolded_string = r"$\bf{BVH}$ $\bf{on}$ $\bf{Intel}$ $\bf{Xeon}$ $\bf{D-2141I}$"
plt.title(bolded_string + "\nRuntime in Cycles", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_cycles.pdf')


plt.close(figure)
plt.plot(triangles, basehovercycles, 'o--', label='Baseline (disjoint)', color='tab:blue')
plt.plot(triangles, basetouchcycles, 'o-', label='Baseline (intersecting)', color='tab:blue')
plt.plot(triangles, vechovercycles, 'x--', label='Vectorized (disjoint)', markersize=10, color='tab:orange')
plt.plot(triangles, vectouchcycles, 'x-', label='Vectorized (intersecting)', markersize=10, color='tab:orange')
plt.plot(triangles, f16chovercycles, '^--', label='Half-Precision (disjoint)', color='tab:green')
plt.plot(triangles, f16ctouchcycles, '^-', label='Half-Precision (intersecting)', color='tab:green')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
a_string = "BVH Runtime in Cycles vs. Input Size "
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Runtime in Cycles vs. Input Size "
bolded_string = r"$\bf{BVH}$ $\bf{on}$ $\bf{Intel}$ $\bf{Xeon}$ $\bf{D-2141I}$"
plt.title(bolded_string + "\nRuntime in Cycles (log. scale)", fontdict=font,loc='left')
plt.yscale('log')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
plt.yticks([50000, 100000, 200000, 500000, 1000000, 2000000, 5000000], ['5×10⁴', '1×10⁵', '2×10⁵', '5×10⁵', '1×10⁶', '2×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_cycles_log.pdf')




basebase = calculateSpeedup(basehovercycles, basehovercycles)
plt.close(figure)
plt.plot(triangles, basebase, 'o-', label='Baseline')
plt.plot(triangles, vechoverspeedup, 'x-', label='Vectorized')
plt.plot(triangles, f16choverspeedup, '^-', label='Half-Precision')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Runtime in Cycles vs. Input Size "
bolded_string = r"$\bf{BVH}$ $\bf{(non-intersecting)}$ $\bf{on}$ $\bf{Intel}$ $\bf{Xeon}$ $\bf{D-2141I}$"
plt.title(bolded_string + "\nSpeedup [Normalized to Baseline]", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_hover_speedup.pdf')

plt.close(figure)
plt.plot(triangles, basebase, 'o-', label='Baseline')
plt.plot(triangles, vectouchspeedup, 'x-', label='Vectorized')
plt.plot(triangles, f16ctouchspeedup, '^-', label='Half-Precision')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Runtime in Cycles vs. Input Size "
bolded_string = r"$\bf{BVH}$ $\bf{(intersecting)}$ $\bf{on}$ $\bf{Intel}$ $\bf{Xeon}$ $\bf{D-2141I}$"
plt.title(bolded_string + "\nSpeedup [Normalized to intersecting or disjoint Baselines]", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_touch_speedup.pdf')




plt.close(figure)
plt.plot(triangles, vechoverspeedup, 'x--', label='Vectorized (disjoint)', color='tab:orange', markersize=10)
plt.plot(triangles, vectouchspeedup, 'x-', label='Vectorized (intersecting)', color='tab:orange', markersize=10)
plt.plot(triangles, f16choverspeedup, '^--', label='Half-Precision (disjoint)', color='tab:green')
plt.plot(triangles, f16ctouchspeedup, '^-', label='Half-Precision (intersecting)', color='tab:green')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Runtime in Cycles vs. Input Size "
bolded_string = r"$\bf{BVH}$ $\bf{Traversal}$ $\bf{on}$ $\bf{Intel}$ $\bf{Xeon}$ $\bf{D-2141I}$"
plt.title(bolded_string + "\nSpeedup [Normalized to intersecting and disjoint Baselines]", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_speedup.pdf')





plt.close(figure)
plt.plot(triangles, basehovertriangletests, 'o--', label='Full Precision (disjoint)', color='deepskyblue')
plt.plot(triangles, basetouchtriangletests, 'o-', label='Full Precision (intersecting)', color='deepskyblue')
plt.plot(triangles, f16chovertriangletests, 'x--', label='Half Precision (disjoint)', color='crimson', markersize=10)
plt.plot(triangles, f16ctouchtriangletests, 'x-', label='Half Precision (intersecting)', color='crimson', markersize=10)
plt.yscale('log')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "Triangle-Triangle Tests performed"
bolded_string = r"$\bf{BVH}$ $\bf{Tree}$ $\bf{Traversal}$"
plt.title(bolded_string + "\nNumber of Triangle-Triangle Tests (log. scale)", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
plt.yticks([10, 100, 1000, 10000, 100000], ['10¹', '10²', '10³', '10⁴', '10⁵'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_triangles.pdf')







plt.close(figure)
plt.plot(triangles, basehoverflops / basehovercycles, 'o--', label='Baseline (disjoint)', color='tab:blue')
plt.plot(triangles, basetouchflops / basetouchcycles, 'o-', label='Baseline (intersecting)', color='tab:blue')
plt.plot(triangles, vechoverflops / vechovercycles, 'x--', label='Vectorized (disjoint)', color='tab:orange', markersize=10)
plt.plot(triangles, vectouchflops / vectouchcycles, 'x-', label='Vectorized (intersecting)', color='tab:orange', markersize=10)
plt.plot(triangles, f16choverflops / f16chovercycles, '^--', label='Half-Precision (disjoint)', color='tab:green')
plt.plot(triangles, f16ctouchflops / f16ctouchcycles, '^-', label='Half-Precision (intersecting)', color='tab:green')
plt.legend()
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size in polygons / 1000', fontdict=font)
plt.grid(axis='y')
a_string = "BVH Performance "
bolded_string = r"$\bf{BVH}$ $\bf{Performance}$ $\bf{Comparison}$"
plt.title(bolded_string + "\nPerformance [flops/cycle]", fontdict=font,loc='left')
plt.xticks(range(0, 5500000, 500000), range(0, 5500, 500), fontsize=16)
#plt.yticks(range(0, 5000001, 1000000), ['0', '1×10⁶', '2×10⁶', '3×10⁶', '4×10⁶', '5×10⁶'], fontsize=16)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.savefig('bvh_performance.pdf')







basehoveropint = basehoverflops / basehovermemory
vechoveropint = vechoverflops / vechovermemory
f16choveropint = f16choverflops / f16chovermemory

basetouchopint = basetouchflops / basetouchmemory
vectouchopint = vectouchflops / vectouchmemory
f16ctouchopint = f16ctouchflops / f16ctouchmemory




basehoverperformance = basehoverflops / basehovercycles
vechoverperformance = vechoverflops / vechovercycles
f16choverperformance = f16choverflops / f16chovercycles

basetouchperformance = basetouchflops / basetouchcycles
vectouchperformance = vectouchflops / vectouchcycles
f16ctouchperformance = f16ctouchflops / f16ctouchcycles




plt.close(figure)
plt.xlim(0.04, 3)
plt.ylim(0.125, 40)
plt.yscale("log", basey=2)
plt.xscale("log", basex=2)

# The memory and compute bound lines for scalar
memoryBound_scalar_x = [0.00404, 0.1294]
memoryBound_scalar_y = [0.125, 4]
computeBound_scalar_x = [0.1294, 8]
computeBound_scalar_y = [4, 4]

# The memory and compute bound lines for AVX2
#memoryBound_avx2_x = [0.1294, 0.5176]
#memoryBound_avx2_y = [4, 16]
#computeBound_avx2_x = [0.5926, 8]
#omputeBound_avx2_y = [16, 16]

# The memory and compute bound lines for AVX512
memoryBound_avx512_x = [0.1294, 1.035]
memoryBound_avx512_y = [4, 32]
computeBound_avx512_x = [1.035, 8]
computeBound_avx512_y = [32, 32]

plt.plot(memoryBound_scalar_x, memoryBound_scalar_y, color='black')
plt.plot(computeBound_scalar_x, computeBound_scalar_y, color='black', label="Theoretical Max. (Scalar)")

#plt.plot(memoryBound_avx2_x, memoryBound_avx2_y, color='orange')
#plt.plot(computeBound_avx2_x, computeBound_avx2_y, "-.", color='orange', label="Theoretical Max. (AVX2)")


plt.plot(memoryBound_avx512_x, memoryBound_avx512_y, color='indigo')
plt.plot(computeBound_avx512_x, computeBound_avx512_y, color='indigo', label="Theoretical Max. (AVX512)")
#plt.plot(instr_mix_avx512_x, instr_mix_avx512_y, "--", color='indigo', label="Instruction Mix (AVX512)")

# Baseline
plt.plot(basehoveropint, basehoverperformance, '.--', color='tab:blue', label='Baseline (disjoint)')  #darkgoldenrod
plt.plot(basetouchopint, basetouchperformance, '.-', color='tab:blue', label='Baseline (intersecting)')

plt.plot(vechoveropint, vechoverperformance, '.--', color='tab:orange', label='Vectorized (disjoint)')  #darkgoldenrod
plt.plot(vectouchopint, vectouchperformance, '.-', color='tab:orange', label='Vectorized (intersecting)')

plt.plot(f16choveropint, f16choverperformance, '.--', color='tab:green', label='Half-Precision (disjoint)')  #darkgoldenrod
plt.plot(f16ctouchopint, f16ctouchperformance, '.-', color='tab:green', label='Half-Precision (intersecting)')

# Optimized Inlined
#plt.plot(gjk_o_i_lu8_op_intenstiy, gjk_o_i_lu8_perf, '.-', color="saddlebrown")

# AoS Vectorized Inlined (AVX2)
#plt.plot(gjk_vectorized_inlined_op_intenstiy, gjk_vectorized_inlined_perf, '.-', color="olive")

# SoA Vectorized (AVX512)
#plt.plot(gjk_soa_vectorized_op_intenstiy, gjk_soa_vectorized_perf, '.-', color="maroon")

#plt.text(0.53,5.45,'$2^{22}$',rotation=0, size=14, color="maroon")
plt.text(1.2,1.1,'Half-Precision',rotation=0, size=12, color="tab:green")
plt.text(0.525,0.9,'Vectorized',rotation=0, size=12, color="tab:orange")
plt.text(0.125,0.5,'Baseline',rotation=0, size=12, color="tab:blue")

#plt.text(0.66,3.2,'$2^{22}$',rotation=0, size=14, color="olive")
#plt.text(1.05,2.87,'AVX2',rotation=0, size=16, color="olive")

#plt.text(0.44,2.13,'$2^{22}$',rotation=0, size=14, color="saddlebrown")
#plt.text(0.6,1.6,'Std. C Opt.',rotation=0, size=16, color="saddlebrown")

#plt.text(0.44,1.3,'$2^{22}$',rotation=0, size=14, color="darkgoldenrod")
#plt.text(0.7,0.98,'Baseline',rotation=0, size=16, color="darkgoldenrod")

#plt.text(3.4,0.63,'$2^4$',rotation=0, size=14, color="maroon")
#plt.text(2.32,22.76,'$2^{13}$',rotation=0, size=14, color="maroon")
#plt.text(2.98,2.94,'$2^{13}$',rotation=0, size=14, color="olive")
#plt.text(1.75,2.05,'$2^{13}$',rotation=0, size=14, color="saddlebrown")
#plt.text(1.28,1.32,'$2^{13}$',rotation=0, size=14, color="darkgoldenrod")

plt.text(0.275,32.75, 'Theoretical Max (AVX512)',rotation=0, size=14, color="indigo")
#plt.text(0.1,16, 'Theoretical Max (AVX2)',rotation=0, size=16, color="darkgoldenrod")
plt.text(0.04,4.15, 'Theoretical Max (Scalar)',rotation=0, size=14, color="black")

# Formatting the plot
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.tick_params(labelsize=16)
plt.xlabel('Operational Intensity [flops/byte]', fontdict=font)
plt.grid(axis='y')
bolded_string = r"$\bf{Roofline}$ $\bf{Plot}$ $\bf{for}$ $\bf{the}$ $\bf{BVH}$ $\bf{Versions}$"
plt.title(bolded_string + '\nPerformance [flops/cycle]', fontdict=font, loc='left')
plt.legend()
# Might need to play with the x value (This line tries to put the legend outside the plot)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.95, top=0.854, wspace=0.2, hspace=0.2)
plt.savefig('bvh_roofline.pdf')















###############################################################################




input_size = [*range(4, 23, 1)]

# Sorry Mihai :(
# It's okay.

# Data points to be used in the plots
#gjk_cycles = [832, 1396, 2451, 4543, 8829, 17058, 33923, 67266, 135177, 271013, 550543, 1083998, 2156256, 4289357, 9093527, 19371503, 39201319, 77078404, 153727260]
#gjk_flops = [942, 1710, 3246, 6318, 12462, 24747, 49323, 98475, 196779, 393390, 786603, 1573035, 3145899, 6291630, 12583083, 22369768, 39146984, 78293936, 134217856]
gjk_cycles = [1288, 1648, 2538, 4646, 8846, 17210, 33884, 67456, 134604, 296952,
              573008, 1092860, 2197684, 4743122, 13658936, 27084038, 38923006, 77632680, 168337672]
gjk_flops = [939, 1707, 3243, 6315, 12459, 24744, 49320, 98472, 196776, 393387,
             786600, 1573032, 3145896, 6291627, 12583080, 25165995, 50331816, 100663464, 201326760]
gjk_perf = calculatePerf(gjk_cycles, gjk_flops)
gjk_speedup = calculateSpeedup(gjk_cycles, gjk_cycles)
gjk_op_intenstiy = calculateOI(gjk_flops)

#gjk_optimized_cycles = [997, 1543, 2710, 5096, 9677, 18944, 37458, 74944, 146550, 296234, 598284, 1172176, 2370467, 4675341, 9843460, 20690487, 41816244, 81566271, 163877846]
#gjk_optimized_flops = [955, 1729, 3264, 6339, 12485, 24772, 49355, 98514, 196821, 393440, 786653, 1573083, 3145951, 6291685, 12583150, 24233944, 46603544, 91342816, 180821264]
#gjk_optimized_cycles = [1518, 1912, 2940, 5112, 9554, 18664, 36814, 73502, 146360, 292468, 583726, 1166012, 2348124, 4937678, 9935122, 20728208, 44535336, 81550414, 163208968]
gjk_optimized_cycles = [1362, 1764, 2688, 4366, 7594, 14196, 27188, 53490,
                        105822, 229104, 420948, 852622, 1751048, 3725272, 7610942, 16095000, 36313424, 61892020, 125989570]
gjk_optimized_flops = [916, 1666, 3225, 6276, 12446, 24709, 49316, 98451,
                       196782, 393377, 786614, 1573020, 3145912, 6291622, 12583111, 25165998, 50331842, 100663463, 201326792]
gjk_optimized_perf = calculatePerf(gjk_optimized_cycles, gjk_optimized_flops)
gjk_optimized_speedup = calculateSpeedup(gjk_optimized_cycles, gjk_cycles)
gjk_optimized_op_intenstiy = calculateOI(gjk_optimized_flops)


#gjk_o_i_cycles = [713, 1103, 1884, 3510, 6496, 12724, 25223, 49433, 102669, 196017, 404012, 789960, 1592501, 3156749, 6907080, 13268662, 25760823, 52499169, 104329888]
#gjk_o_i_flops = [955, 1729, 3264, 6339, 12485, 24772, 49355, 98514, 196821, 393440, 786653, 1573083, 3145951, 6291685, 12583150, 24233944, 46603544, 91342816, 180821264]
gjk_o_i_cycles = [1386, 1530, 2222, 3650, 6770, 12820, 25108, 49494,
                  99038, 196560, 404182, 789498, 1577400, 3245516, 6545216, 13216116, 31206406, 52733144, 104070992]
gjk_o_i_flops = [939, 1707, 3243, 6315, 12459, 24744, 49320, 98472,
                 196776, 393387, 786600, 1573032, 3145896, 6291627, 12583080, 25165995, 50331816, 100663464, 201326760]
gjk_o_i_perf = calculatePerf(gjk_o_i_cycles, gjk_o_i_flops)
gjk_o_i_speedup = calculateSpeedup(gjk_o_i_cycles, gjk_cycles)
gjk_o_i_op_intenstiy = calculateOI(gjk_o_i_flops)

#gjk_o_i_lu8_cycles = [703, 1078, 1925, 3537, 6600, 12769, 25430, 49976, 102125, 201543, 401974, 784403, 1672933, 3717257, 7245230, 13559705, 25909261, 56818646, 103202148]
#gjk_o_i_lu8_flops = [957, 1731, 3264, 6338, 12485, 24771, 49353, 98513, 196820, 393442, 786656, 1573088, 3145957, 6291691, 12583159, 23767920, 44739432, 86682464, 163298400]
#gjk_o_i_lu8_cycles = [1034, 1364, 2070, 3594, 6594, 12690, 24744, 49182, 98492, 204464, 392926, 790148, 1583186, 3236690, 6588126, 13570806, 30663150, 51680260, 103158384]
gjk_o_i_lu8_cycles = [1244, 1442, 2168, 3522, 6226, 11816, 22788, 44804, 89442,
                      178064, 355698, 711398, 1438122, 2966340, 6119580, 12547814, 29151156, 47909034, 95056348]
gjk_o_i_lu8_flops = [939, 1707, 3243, 6315, 12459, 24744, 49320, 98472, 196776,
                     393387, 786600, 1573032, 3145896, 6291627, 12583080, 25165995, 50331816, 100663464, 201326760]
gjk_o_i_lu8_perf = calculatePerf(gjk_o_i_lu8_cycles, gjk_o_i_lu8_flops)
gjk_o_i_lu8_speedup = calculateSpeedup(gjk_o_i_lu8_cycles, gjk_cycles)
gjk_o_i_lu8_op_intenstiy = calculateOI(gjk_o_i_lu8_flops)

#gjk_vectorized_cycles = [1032, 1185, 1813, 3183, 5947, 11565, 22694, 45438, 89778, 180614, 358988, 717051, 1440319, 2884327, 6612599, 13267698, 25486536, 51851147, 102994038]
#gjk_vectorized_flops = [1454, 2606, 4910, 9518, 18734, 37163, 74027, 147755, 295211, 590126, 1179947, 2359595, 4718891, 9437486, 18874652, 37748984, 76429792, 167857136, 268435456]
gjk_vectorized_cycles = [1308, 1538, 2302, 3942, 6524, 12110, 23190, 45120, 89398,
                         212304, 385456, 780596, 1544352, 3066462, 6371484, 13222392, 30793866, 52064938, 103457734]
gjk_vectorized_flops = [1323, 2347, 4395, 8491, 16683, 33064, 65832, 131368, 262440,
                        524587, 1048872, 2097448, 4194600, 8388907, 16777512, 33554731, 67109160, 134218024, 268435752]
gjk_vectorized_perf = calculatePerf(gjk_vectorized_cycles, gjk_vectorized_flops)
gjk_vectorized_speedup = calculateSpeedup(gjk_vectorized_cycles, gjk_cycles)
gjk_vectorized_op_intenstiy = calculateOI(gjk_vectorized_flops)

#gjk_vectorized_inlined_cycles = [719, 1014, 1693, 3066, 5802, 11291, 22224, 44558, 88453, 176773, 356066, 706017, 1417855, 2858754, 5903309, 11889467, 23199334, 46924311, 94008049]
#gjk_vectorized_inlined_flops = [1390, 2478, 4654, 9006, 17710, 35115, 69931, 139563, 278827, 557358, 1114411, 2228523, 4456747, 8913198, 17826084, 35651864, 72043664, 155929728, 268435456]
gjk_vectorized_inlined_cycles = [1118, 1354, 2262, 3712, 6582, 12140, 23562, 45940, 92924,
                                 217814, 390120, 768386, 1500494, 3019438, 6298546, 12161234, 26076978, 48472878, 96848942]
gjk_vectorized_inlined_flops = [1323, 2347, 4395, 8491, 16683, 33064, 65832, 131368, 262440,
                                524587, 1048872, 2097448, 4194600, 8388907, 16777512, 33554731, 67109160, 134218024, 268435752]
gjk_vectorized_inlined_perf = calculatePerf(gjk_vectorized_inlined_cycles, gjk_vectorized_inlined_flops)
gjk_vectorized_inlined_speedup = calculateSpeedup(gjk_vectorized_inlined_cycles, gjk_cycles)
gjk_vectorized_inlined_op_intenstiy = calculateOI(gjk_vectorized_inlined_flops)

##############################
#gjk_soa_v_slow_idx_cycles = [952, 1081, 1500, 2583, 4609, 8732, 16987, 33546, 66724, 132354, 264412, 531740, 1055783, 2125153, 4465748, 8979228, 18490471, 36086785, 72229488]
#gjk_soa_v_slow_idx_flops = [1326, 2350, 4398, 8494, 16686, 33067, 65835, 131371, 262443, 524590, 1048875, 2097451, 4194603, 8388910, 16777508, 33554680, 67109104, 134217968, 268435584]
gjk_soa_v_slow_idx_cycles = [1364, 1250, 1660, 2806, 4954, 9012, 17160, 33442, 66660,
                             132762, 295854, 574254, 1193932, 2281082, 4712666, 11130410, 20698188, 43142480, 72427146]
gjk_soa_v_slow_idx_flops = [1195, 2091, 3883, 7467, 14635, 28968, 57640, 114984, 229672,
                            459051, 917800, 1835304, 3670312, 7340331, 14680360, 29360427, 58720552, 117440808, 234881320]
gjk_soa_v_slow_idx_speedup = calculateSpeedup(gjk_soa_v_slow_idx_cycles, gjk_cycles)
gjk_soa_v_slow_idx_perf = calculatePerf(gjk_soa_v_slow_idx_cycles, gjk_soa_v_slow_idx_flops)

#gjk_soa_v_slow_idx_inl_cycles = [741, 842, 1048, 1619, 2720, 4905, 9255, 18072, 35788, 71062, 141867, 282679, 565195, 1124478, 2382808, 5203033, 10719750, 22378936, 45572476]
#gjk_soa_v_slow_idx_inl_flops = [1326, 2350, 4398, 8494, 16686, 33067, 65835, 131371, 262443, 524590, 1048875, 2097451, 4194603, 8388910, 16777508, 33554680, 67109104, 134217968, 268435584]
gjk_soa_v_slow_idx_inl_cycles = [1682, 1154, 1480, 1900, 2932, 5034, 9360, 18092, 36264,
                                 71782, 174804, 327340, 657634, 1322302, 2734544, 6169992, 12855090, 24414858, 45948272]
gjk_soa_v_slow_idx_inl_flops = [1195, 2091, 3883, 7467, 14635, 28968, 57640, 114984, 229672,
                                459051, 917800, 1835304, 3670312, 7340331, 14680360, 29360427, 58720552, 117440808, 234881320]
gjk_soa_v_slow_idx_inl_speedup = calculateSpeedup(gjk_soa_v_slow_idx_inl_cycles, gjk_cycles)
gjk_soa_v_slow_idx_inl_perf = calculatePerf(gjk_soa_v_slow_idx_inl_cycles, gjk_soa_v_slow_idx_inl_flops)
###############################

#gjk_soa_vectorized_cycles = [917, 1026, 1084, 1301, 1591, 2265, 3510, 6179, 11450, 22061, 84928, 157155, 283058, 668896, 1500451, 5234456, 11749913, 25483926, 52995337]
#gjk_soa_vectorized_flops = [1326, 2350, 4398, 8494, 16686, 33067, 65835, 131371, 262443, 524590, 1048875, 2097451, 4194603, 8388910, 16777508, 33554680, 67109104, 134217968, 268435584]
gjk_soa_vectorized_cycles = [1610, 1334, 1340, 1614, 1752, 2644, 3576, 6114, 12220,
                             22618, 121818, 213858, 425502, 914462, 2034442, 6162966, 19821088, 27845568, 54005966]
gjk_soa_vectorized_flops = [1195, 2091, 3883, 7467, 14635, 28968, 57640, 114984, 229672,
                            459051, 917800, 1835304, 3670312, 7340331, 14680360, 29360427, 58720552, 117440808, 234881320]
gjk_soa_vectorized_perf = calculatePerf(gjk_soa_vectorized_cycles, gjk_soa_vectorized_flops)
gjk_soa_vectorized_speedup = calculateSpeedup(gjk_soa_vectorized_cycles, gjk_cycles)
gjk_soa_vectorized_op_intenstiy = calculateOI(gjk_soa_vectorized_flops)

#gjk_soa_vectorized_lu2_cycles = [868, 1039, 1123, 1312, 1628, 2271, 3560, 6127, 11394, 22087, 75752, 152408, 331819, 666783, 1676025, 4962694, 12361912, 25375337, 53995497]
#gjk_soa_vectorized_lu2_flops = [942, 2606, 4654, 8750, 16942, 33323, 66091, 131627, 262699, 524846, 1049131, 2097707, 4194859, 8389166, 16777748, 33554928, 67109344, 134218208, 268435808]
gjk_soa_vectorized_lu2_cycles = [1270, 1840, 1206, 1496, 1846, 2820, 3602, 6098, 12606,
                                 23308, 113112, 181706, 407890, 966388, 2735116, 6455846, 20262788, 27752144, 54090402]
gjk_soa_vectorized_lu2_flops = [939, 2347, 4139, 7723, 14891, 29224, 57896, 115240, 229928,
                                459307, 918056, 1835560, 3670568, 7340587, 14680616, 29360683, 58720808, 117441064, 234881576]
gjk_soa_vectorized_lu2_perf = calculatePerf(gjk_soa_vectorized_lu2_cycles, gjk_soa_vectorized_lu2_flops)
gjk_soa_vectorized_lu2_speedup = calculateSpeedup(gjk_soa_vectorized_lu2_cycles, gjk_cycles)
gjk_soa_vectorized_lu2_op_intenstiy = calculateOI(gjk_soa_vectorized_lu2_flops)

#gjk_soa_vectorized_inlined_cycles = [740, 819, 866, 1074, 1318, 1898, 3112, 5410, 10822, 21335, 79923, 148984, 305713, 615233, 1603444, 4306188, 10642093, 23584374, 50231125]
#gjk_soa_vectorized_inlined_flops = [1262, 2222, 4142, 7982, 15662, 31019, 61739, 123179, 246059, 491822, 983339, 1966379, 3932459, 7864622, 15728939, 31457568, 62914840, 125829424, 251658528]
gjk_soa_vectorized_inlined_cycles = [1704, 1088, 964, 1108, 1586, 1970, 3124, 5660, 12112,
                                     21920, 119996, 203478, 384966, 907782, 1955504, 5476374, 18170646, 25809308, 49872214]
gjk_soa_vectorized_inlined_flops = [1195, 2091, 3883, 7467, 14635, 28968, 57640, 114984,
                                    229672, 459051, 917800, 1835304, 3670312, 7340331, 14680360, 29360427,
                                    58720552, 117440808, 234881320]
gjk_soa_vectorized_inlined_perf = calculatePerf(gjk_soa_vectorized_inlined_cycles, gjk_soa_vectorized_inlined_flops)
gjk_soa_vectorized_inlined_speedup = calculateSpeedup(gjk_soa_vectorized_inlined_cycles, gjk_cycles)
gjk_soa_vectorized_inlined_op_intenstiy = calculateOI(gjk_soa_vectorized_inlined_flops)

###################################### Speedup Plots #########################################################################
'''
plt.plot(input_size, gjk_speedup, 'o-', label='Baseline')
plt.plot(input_size, gjk_optimized_speedup, 'o-', label='Optimized')
#plt.plot(input_size, gjk_o_i_speedup, 'o-', label='Optimized Inlined')
plt.plot(input_size, gjk_o_i_lu8_speedup, 'o-', label='O. I. Loop Unroll (LU) 8')
#plt.plot(input_size, gjk_vectorized_speedup, 'o-', label='AoS Vectorized (AVX2)')
plt.plot(input_size, gjk_vectorized_inlined_speedup, 'o-', label='AoS Vectorized Inlined (AVX2)')
#plt.plot(input_size, gjk_soa_v_slow_idx_speedup, 'o-', label='SoA Vectorized Slow Idx Calc. (AVX512)')
plt.plot(input_size, gjk_soa_v_slow_idx_inl_speedup, 'o-', label='SoA Vectorized S.I.C Inl. (AVX512)')
#plt.plot(input_size, gjk_soa_vectorized_speedup, 'o-', label='SoA Vectorized (AVX512)')
#plt.plot(input_size, gjk_soa_vectorized_lu2_speedup, 'o-', label='SoA Vectorized LU 2 (AVX512)')
plt.plot(input_size, gjk_soa_vectorized_inlined_speedup, 'o-', label='SoA Vectorized Inlined (AVX512)')

font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.xlabel('Input Size [n] (n is the power in $2^n$)', fontdict=font)
plt.grid(axis='y')
a_string = "GJK Speedup vs Input Size"
# Make string bold
bolded_string = r"$\bf{GJK}$ $\bf{(single}$ $\bf{precision)}$ $\bf{on}$ $\bf{Intel}$ $\bf{Core}$ $\bf{i7}$"
plt.title(bolded_string + "\nSpeedup [Normalized to Baseline]", fontdict=font,loc='left')

# Labels on the plot
plt.text(11.5,13,'L1',rotation=0, size=16)
plt.text(15,13,'L2',rotation=0, size=16)
plt.text(19,13,'L3',rotation=0, size=16)
plt.text(10.25,14.36,'SoA Strength Reduction 13.5x',rotation=0, size=16, color='brown')
plt.text(8 ,4.7,'SoA AVX512 4.1x',rotation=0, size=16, color='indigo')
plt.text(8.45, 2,'Std C Inlining / AoS AVX2 1.5x',rotation=0, size=16, color='maroon')
plt.text(8.45, 0.3,'Baseline 1x / Block Optimizations 1.25x',rotation=0, size=16, color='darkgoldenrod')

plt.ylim(0, 15)
plt.yticks(range(0, 16, 2), fontsize=16)
plt.xticks(range(4, 23, 1), fontsize=16)

# Cache lines
plt.vlines(x=11, ymin=0, ymax=15, colors='gray')
plt.vlines(x=14.415, ymin=0, ymax=15, colors='gray')
plt.vlines(x=18.415, ymin=0, ymax=15, colors='gray')

figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.95, top=0.854, wspace=0.2, hspace=0.2)
plt.savefig('gjk_speedup.pdf')
#plt.show()
'''
###################################### Performance Plots (Standard C) #####################################################################
'''
plt.plot(input_size, gjk_perf, '.-', color="darkgoldenrod")
plt.plot(input_size, gjk_optimized_perf, '.-', color="saddlebrown")
plt.plot(input_size, gjk_o_i_perf, '.-', color="olive")

font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}

plt.xlabel('Input Size [n] (n is the power in $2^n$)', fontdict=font)
plt.grid(axis='y')
bolded_string = r"$\bf{Baseline}$ $\bf{/}$ $\bf{Standard}$ $\bf{C}$ $\bf{Optimizations}$ $\bf{Performance}$ $\bf{Comparison}$"
plt.title(bolded_string + '\nPerformance [flops/cycle]', fontdict=font,loc='left')
plt.ylim(0, 4)
plt.xlim(3, 23)
plt.yticks(range(0, 5, 1), size=16)
plt.xticks(range(4, 23, 1), size=16)

# Labels on the plot
plt.text(11.5,3.2,'L1',rotation=0, size=16)
plt.text(15,3.2,'L2',rotation=0, size=16)
plt.text(19,3.2,'L3',rotation=0, size=16)
plt.text(4, 2.2,'Std C Inlining / Joint Loop Optimization',rotation=0, size=16, color="olive")
plt.text(8.6, 1.10,'Baseline',rotation=0, size=16, color="darkgoldenrod")
plt.text(8.6, 1.61,'Block Optimizations',rotation=0, size=16, color="saddlebrown")
plt.text(4, 2.6,'Max. Perf. (Dependencies)',rotation=0, size=16)
plt.text(4, 3.1,'Max. Perf. (Instruction Mix)',rotation=0, size=16)

# Peak performance lines
plt.hlines(xmin=3, xmax=23, y=3, linestyles='dashed', colors='gray')
plt.hlines(xmin=3, xmax=23, y=2.52, linestyles='dashed', colors='gray')

# Cache lines
plt.vlines(x=11, ymin=0, ymax=4, colors='gray')
plt.vlines(x=14.415, ymin=0, ymax=4, colors='gray')
plt.vlines(x=18.415, ymin=0, ymax=4, colors='gray')

figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.95, top=0.854, wspace=0.2, hspace=0.2)
plt.savefig('gjk_performance_stdc.pdf')
'''
###################################### Performance Plots (Vectorized) #####################################################################

'''
plt.plot(input_size, gjk_vectorized_inlined_perf, '.-', color="darkgoldenrod")
plt.plot(input_size, gjk_soa_v_slow_idx_inl_perf, '.-', color="saddlebrown")
plt.plot(input_size, gjk_soa_vectorized_perf, '.-', color="olive")
plt.plot(input_size, gjk_soa_vectorized_inlined_perf, '.-', color="maroon")

font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}

plt.xlabel('Input Size [n] (n is the power in $2^n$)', fontdict=font)
plt.grid(axis='y')
bolded_string = r"$\bf{AVX2}$ $\bf{/}$ $\bf{AVX512}$ $\bf{Performance}$ $\bf{Comparison}$"
plt.title(bolded_string + '\nPerformance [flops/cycle]', fontdict=font,loc='left')
plt.ylim(0, 32)
plt.xlim(3, 23)
plt.yticks(range(0, 33, 2), fontsize=16)
plt.xticks(range(4, 23, 1), fontsize=16)

plt.text(11.5,28.5,'L1',rotation=0, size=16)
plt.text(15,28.5,'L2',rotation=0, size=16)
plt.text(19,28.5,'L3',rotation=0, size=16)
plt.text(8.6, 4.37, 'SoA V. Inl. (AVX512)',rotation=0, size=16, color="saddlebrown")
plt.text(8.6, 1,'AoS Vectorized Inlined (AVX2)',rotation=0, size=16, color="darkgoldenrod")
plt.text(7.72, 6.64,'SoA V. Strength Reduction',rotation=0, size=16, color="olive")
plt.text(5.7, 21,'SoA V. S.R. Inl. (AVX512)',rotation=0, size=16, color="maroon")
plt.text(14.5, 23,'Max. Perf. (Instruction Mix)',rotation=0, size=16)

plt.vlines(x=11, ymin=0, ymax=32, colors='gray')
plt.vlines(x=14.415, ymin=0, ymax=32, colors='gray')
plt.vlines(x=18.415, ymin=0, ymax=32, colors='gray')
plt.hlines(xmin=3, xmax=23, y=22.4, linestyles='dashed', colors='gray')

figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.95, top=0.854, wspace=0.2, hspace=0.2)
plt.savefig('gjk_performance_avx.pdf')
'''
###################################### Roofline Plot #####################################################################
'''
plt.xlim(0.01923, 8)
plt.ylim(0.5, 40)
plt.yscale("log", base=2)
plt.xscale("log", base=2)

# The memory and compute bound lines for scalar
memoryBound_scalar_x = [0.0185, 0.148]
memoryBound_scalar_y = [0.5, 4]
computeBound_scalar_x = [0.148, 8]
computeBound_scalar_y = [4, 4]

# The memory and compute bound lines for AVX2
memoryBound_avx2_x = [0.148, 0.5926]
memoryBound_avx2_y = [4, 16]
computeBound_avx2_x = [0.5926, 8]
computeBound_avx2_y = [16, 16]

# The memory and compute bound lines for AVX512
memoryBound_avx512_x = [0.5926, 1.185]
memoryBound_avx512_y = [16, 32]
computeBound_avx512_x = [1.185, 8]
computeBound_avx512_y = [32, 32]

plt.plot(memoryBound_scalar_x, memoryBound_scalar_y, color='black')
plt.plot(computeBound_scalar_x, computeBound_scalar_y, color='black', label="Theoretical Max. (Scalar)")

plt.plot(memoryBound_avx2_x, memoryBound_avx2_y, color='orange')
plt.plot(computeBound_avx2_x, computeBound_avx2_y, "-.", color='orange', label="Theoretical Max. (AVX2)")


plt.plot(memoryBound_avx512_x, memoryBound_avx512_y, color='indigo')
plt.plot(computeBound_avx512_x, computeBound_avx512_y, color='indigo', label="Theoretical Max. (AVX512)")
#plt.plot(instr_mix_avx512_x, instr_mix_avx512_y, "--", color='indigo', label="Instruction Mix (AVX512)")

# Baseline
plt.plot(gjk_op_intenstiy, gjk_perf, '.-', color="darkgoldenrod")

# Optimized Inlined
plt.plot(gjk_o_i_lu8_op_intenstiy, gjk_o_i_lu8_perf, '.-', color="saddlebrown")

# AoS Vectorized Inlined (AVX2)
plt.plot(gjk_vectorized_inlined_op_intenstiy, gjk_vectorized_inlined_perf, '.-', color="olive")

# SoA Vectorized (AVX512)
plt.plot(gjk_soa_vectorized_op_intenstiy, gjk_soa_vectorized_perf, '.-', color="maroon")

plt.text(0.53,5.45,'$2^{22}$',rotation=0, size=14, color="maroon")
plt.text(1.27,10.84,'AVX512',rotation=0, size=16, color="maroon")

plt.text(0.66,3.2,'$2^{22}$',rotation=0, size=14, color="olive")
plt.text(1.05,2.87,'AVX2',rotation=0, size=16, color="olive")

plt.text(0.44,2.13,'$2^{22}$',rotation=0, size=14, color="saddlebrown")
plt.text(0.6,1.6,'Std. C Opt.',rotation=0, size=16, color="saddlebrown")

plt.text(0.44,1.3,'$2^{22}$',rotation=0, size=14, color="darkgoldenrod")
plt.text(0.7,0.98,'Baseline',rotation=0, size=16, color="darkgoldenrod")

plt.text(3.4,0.63,'$2^4$',rotation=0, size=14, color="maroon")
plt.text(2.32,22.76,'$2^{13}$',rotation=0, size=14, color="maroon")
plt.text(2.98,2.94,'$2^{13}$',rotation=0, size=14, color="olive")
plt.text(1.75,2.05,'$2^{13}$',rotation=0, size=14, color="saddlebrown")
plt.text(1.28,1.32,'$2^{13}$',rotation=0, size=14, color="darkgoldenrod")

plt.text(0.96,33.8, 'Theoretical Max (AVX512)',rotation=0, size=16, color="indigo")
plt.text(0.1,16, 'Theoretical Max (AVX2)',rotation=0, size=16, color="darkgoldenrod")
plt.text(0.02,4, 'Theoretical Max (Scalar)',rotation=0, size=16, color="black")

# Formatting the plot
font = {
    'family': 'sans serif',
    'color': 'black',
    'weight': 'normal',
    'size': 18
}
plt.tick_params(labelsize=16)
plt.xlabel('Operational Intensity [flops/byte]', fontdict=font)
plt.grid(axis='y')
bolded_string = r"$\bf{Roofline}$ $\bf{Plot}$ $\bf{for}$ $\bf{Various}$ $\bf{GJK}$ $\bf{Versions}$"
plt.title(bolded_string + '\nPerformance [flops/cycle]', fontdict=font, loc='left')

# Might need to play with the x value (This line tries to put the legend outside the plot)
figure = plt.gcf()
figure.set_size_inches(10, 7)
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.95, top=0.854, wspace=0.2, hspace=0.2)
plt.savefig('gjk_roofline.pdf')
'''
