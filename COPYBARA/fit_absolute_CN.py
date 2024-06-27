from multiprocessing import Pool
from multiprocessing import cpu_count
import numpy as np
from scipy.stats import gaussian_kde
# from scipy.stats import skew as scipy_skew
# from scipy.stats import kurtosis as scipy_kurtosis

import CN_helper as cnh

#----
# 1. Parse command line arguments
#----


#----
# 2. Define functions
#----
# def is_constant(x):
#     return np.all(x == x[0])

# def Modes(x, min_size=0.1):
#     '''
#     The Modes function is a simple, deterministic function that differences the kernel density of x and reports a number of modes equal to half the number of changes in direction, 
#     although the min.size function can be used to reduce the number of modes returned, and defaults to 0.1, eliminating modes that do not have at least 10% of the distributional area. 
#     The Modes function returns a list with three components: modes, modes.dens, and size. 
#     The elements in each component are ordered according to the decreasing density of the modes. The modes component is a vector of the values of x associated with the modes. 
#     The modes.dens component is a vector of the kernel density estimates at the modes. 
#     The size component is a vector of the proportion of area underneath each mode. 
#     Function adapted from the Modes() function of the "LaplacesDemon" R package
#     '''
#     if x is None:
#         raise ValueError("The x argument is required.")
#     # Convert input to a numpy array of floats and filter out non-finite values
#     x = np.asarray(x).astype(float)
#     x = x[np.isfinite(x)]
#     # Check if the input array is constant
#     if is_constant(x):
#         return {'modes': [np.nan], 'mode_dens': [np.nan], 'size': [1]}
#     # Compute the kernel density estimate
#     density = gaussian_kde(x)
#     # Create a grid of 1000 points over the range of the data
#     x_grid = np.linspace(np.min(x), np.max(x), 1000)
#     # Evaluate the density over the grid
#     dens_y = density(x_grid)
#     # Compute the differences between consecutive density values
#     dens_y_diff = np.diff(dens_y)
#     # Identify where the density is increasing (1) or not increasing (0)
#     incr = np.where(dens_y_diff > 0, 1, 0)
#     # Initialize the list of segment boundaries
#     begin = [0]
#     for i in range(1, len(incr)):
#         if incr[i] != incr[i - 1]:
#             begin.append(i)
#     # Add the end of the array as the final boundary
#     begin.append(len(incr))
#     #
#     size = []
#     modes = []
#     mode_dens = []
#     # Sum of all density values for normalization
#     dens_y_sum = np.sum(dens_y)
#     #
#     j = 0
#     while j < len(begin) - 1:
#         start_idx = begin[j] # Define the start and end of the current segment
#         end_idx = begin[j + 2] if j + 2 < len(begin) else len(dens_y) # Extract the segment of the density
#         segment = dens_y[start_idx:end_idx]  # Calculate the size of the segment
#         segment_size = np.sum(segment) / dens_y_sum
#         if segment_size >= min_size: # Define the grid points and density values for the segment
#             kde_x = x_grid[start_idx:end_idx]
#             kde_y = segment
#             mode_idx = np.argmax(kde_y) # Find the index of the mode (maximum density) within the segment
#             modes.append(kde_x[mode_idx]) # Store the mode and its density
#             mode_dens.append(kde_y[mode_idx])
#             size.append(segment_size) # Store the size of the segment
#         j += 2
#     # Order the results by density in descending order
#     order = np.argsort(mode_dens)[::-1]
#     size = np.array(size)[order]
#     modes = np.array(modes)[order]
#     mode_dens = np.array(mode_dens)[order]
#     # Normalize the sizes to sum to 1 if their sum exceeds 1
#     if np.sum(size) > 1:
#         size = size / np.sum(size)
#     # Return dictionary
#     return {'modes': modes, 'mode_dens': mode_dens, 'size': size}

# def is_unimodal(x, min_size=0.1):
#     '''
#     Function to test if data is unimodal
#     '''
#     modes_result = Modes(x, min_size)
#     if len(modes_result['modes']) == 1 and not np.any(np.isnan(modes_result['modes'])):
#         return True
#     else:
#         return False

# def skew(x, na_rm=True, type=3):
#     '''
#     Calculate the skewness of the input array x based on the specified type.
#     '''
#     if na_rm:
#         x = x[~np.isnan(x)]
#     n = len(x)
#     if n == 0:
#         return np.nan
#     #
#     if type == 1:
#         skewer = np.sqrt(n) * (np.sum((x - np.mean(x))**3) / (np.sum((x - np.mean(x))**2)**(3/2)))
#     elif type == 2:
#         skewer = n * np.sqrt(n - 1) * (np.sum((x - np.mean(x))**3) / ((n - 2) * np.sum((x - np.mean(x))**2)**(3/2)))
#     elif type == 3:
#         skewer = np.sum((x - np.mean(x))**3) / (n * np.std(x)**3)
#     else:
#         raise ValueError("Type must be 1, 2, or 3")
#     #
#     return skewer

# def kurtosi(x, na_rm=True, type=3):
#     '''
#     Calculate the kurtosis of the input array x based on the specified type.
#     '''
#     if na_rm:
#         x = x[~np.isnan(x)]
#     n = len(x)
#     if n == 0:
#         return np.nan
#     #
#     if type == 1:
#         kurt = (np.sum((x - np.mean(x))**4) * n) / (np.sum((x - np.mean(x))**2)**2) - 3
#     elif type == 2:
#         kurt = (n * (n + 1) * np.sum((x - np.mean(x))**4)) / ((n - 1) * (n - 2) * (n - 3) * (np.sum((x - np.mean(x))**2) / (n - 1))**2) - 3 * (n - 1)**2 / ((n - 2) * (n - 3))
#     elif type == 3:
#         kurt = (np.sum((x - np.mean(x))**4) / (n * np.std(x)**4)) - 3
#     else:
#         raise ValueError("Type must be 1, 2, or 3")
#     #
#     return kurt


# def bimodality_coefficient(x, na_rm=False):
#     '''
#     Calculate the bimodality coefficient using the skew and kurtosi functions.
#     '''
#     if na_rm:
#         x = x[~np.isnan(x)]
#     n = len(x)
#     if n == 0:
#         return np.nan
#     #
#     m3 = skew(x, na_rm=na_rm, type=2)
#     # m3 = scipy_skew(x)
#     m4 = kurtosi(x, na_rm=na_rm, type=2)
#     # m4 = scipy_kurtosis(x)
#     bc = (m3**2 + 1) / (m4 + 3 * ((n - 1)**2 / ((n - 2) * (n - 3))))
#     return bc


# def mode(x):
#     if x is None: 
#         print("The x argument in Modes function is required.")
#         return None
#     npx = np.asarray(x)[np.isfinite(np.asarray(x))]
#     if is_constant(npx):
#         return np.nan
#     if np.all(npx == np.round(npx)):
#         values, counts = np.unique(npx, return_counts=True)
#         mode_value = values[np.argmax(counts)]
#     else:
#         kde = gaussian_kde(npx)
#         kde_values = kde.evaluate(npx)
#         mode_value = x[np.argmax(kde_values)]
#     return mode_value


#----
# 3. Estimate purity using phased heterozygous SNPs
#----
# Method/principles based on: https://doi.org/10.1371/journal.pone.0045835
allele_counts_path = "chr17_215000645_allele_counts_hetSNPs.tsv"
# process input - allele_counts_hetSNPs.tsv
allele_counts = []
with open(allele_counts_path, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        A,C,G,T,N = int(fields[5]),int(fields[6]),int(fields[7]),int(fields[8]),int(fields[9])
        DP = A+C+G+T+N
        # only include het SNPs allele counts with a total depth of > 20.
        if DP > 20:
            fields.append(DP)
            # print(fields.append(str(DP)))
            allele_counts.append(fields)

# columns:"chr", "start", "end","REF", "ALT", "A","C","G","T","N","AF_0","AF_1","GT","PS","DP"

# get unique phase sets to iterate (or multiprocess through)
phasesets = list(dict.fromkeys([x[-2] for x in allele_counts]))
ps = '72909885'
ac_ps = [x for x in allele_counts if x[-2] == ps]




# filter snps with DP<20
# get mean DP for each phasing set
# remove SNPs with AF_0 or AF_1 = 1 or = 0


# data = [1, 2, 2, 3, 3, 3, 4, 4, 10, 20, 20, 30, 30, 30, 40, 40]
# test = [0.3068182, 0.3522727, 0.2873563, 0.2926829, 0.4230769, 0.2676056, 0.3, 0.15625, 0.2459016, 0.3703704, 0.3157895, 0.2258065, 0.2419355, 0.265625, 0.2686567, 0.234375, 0.2537313, 0.2461538, 0.1343284, 0.3015873, 0.2205882, 0.2173913, 0.2878788, 0.203125, 0.2537313, 0.1940299, 0.2461538, 0.2083333, 0.2191781, 0.2266667, 0.2222222, 0.2575758, 0.1666667, 0.2368421, 0.1025641, 0.2463768, 0.2266667, 0.3670886, 0.225, 0.4, 0.1627907, 0.2289157, 0.3085106, 0.2857143, 0.2926829, 0.2, 0.278481, 0.4727273, 0.3214286, 0.2183908, 0.2592593, 0.3544304, 0.3291139, 0.5633803, 0.4025974, 0.3026316, 0.2962963, 0.2191781, 0.3066667, 0.304878, 0.2771084, 0.2597403, 0.3658537, 0.3076923, 0.2643678, 0.3295455, 0.3555556, 0.3222222, 0.3095238, 0.3181818, 0.3209877, 0.28125, 0.2717391, 0.2527473, 0.3870968, 0.2666667, 0.25, 0.255814, 0.2674419, 0.2467532, 0.297619, 0.3037975, 0.3493976, 0.443038, 0.3846154, 0.4558824, 0.375, 0.28, 0.2739726, 0.3, 0.3478261, 0.3484848, 0.3428571, 0.40625, 0.2786885, 0.6896552, 0.3823529, 0.3043478, 0.28125, 0.4067797, 0.6590909, 0.6363636, 0.7126437, 0.6585366, 0.5769231, 0.7183099, 0.7, 0.828125, 0.7540984, 0.6111111, 0.6666667, 0.7741935, 0.7258065, 0.71875, 0.7313433, 0.6875, 0.7462687, 0.6461538, 0.8507463, 0.6984127, 0.7647059, 0.7391304, 0.6818182, 0.796875, 0.7164179, 0.7910448, 0.7230769, 0.7777778, 0.739726, 0.7066667, 0.7638889, 0.7272727, 0.8333333, 0.7368421, 0.8717949, 0.7246377, 0.7733333, 0.6202532, 0.75, 0.5882353, 0.8255814, 0.7710843, 0.6914894, 0.7142857, 0.695122, 0.7866667, 0.7088608, 0.5272727, 0.6547619, 0.7816092, 0.7283951, 0.6455696, 0.6455696, 0.4366197, 0.5844156, 0.6973684, 0.691358, 0.7671233, 0.68, 0.6585366, 0.686747, 0.7012987, 0.6219512, 0.6794872, 0.7241379, 0.6590909, 0.6222222, 0.6666667, 0.6904762, 0.6590909, 0.6666667, 0.7083333, 0.7065217, 0.7032967, 0.5913978, 0.6666667, 0.7291667, 0.7325581, 0.7209302, 0.7402597, 0.702381, 0.6835443, 0.6506024, 0.556962, 0.6025641, 0.5441176, 0.5972222, 0.7066667, 0.7260274, 0.6714286, 0.6521739, 0.6363636, 0.6, 0.59375, 0.6557377, 0.3103448, 0.5735294, 0.6521739, 0.671875, 0.5932203]

# print(cnh.is_unimodal(test))
# print(cnh.bimodality_coefficient(test))


# print(cnh.is_unimodal(data))
# print(cnh.bimodality_coefficient(data))