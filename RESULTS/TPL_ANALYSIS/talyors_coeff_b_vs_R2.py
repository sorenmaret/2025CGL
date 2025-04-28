from scipy.stats import linregress, nbinom, lognorm
import numpy as np
import matplotlib.pyplot as plt

# Define a small epsilon to prevent log10 of zero or negative values
EPSILON = 1e-10

# Define the means and group size

p = 0.3
n_values = np.linspace(10, 1000, 30, dtype=int)
b_values_binom, r2_values_binom = [], []

means = np.logspace(0, 2, 10)  # logarithmically spaced means from 1 to 100
group_size = 1000  # number of samples per group

# We'll vary dispersion (r) for Negative Binomial and shape (sigma) for Log-Normal
dispersion_values = np.linspace(1, 50, 30)  # for NB: smaller = more overdispersed
sigma_values = np.linspace(0.2, 2.5, 30)    # for Log-Normal: larger = more skewed

b_values_nb, r2_values_nb = [], []
b_values_ln, r2_values_ln = [], []

b_values_gauss, r2_values_gauss = [], []
std_devs = np.linspace(1, 20, 30)

# Simulate Taylor's law for each value of dispersion/skew

for n in n_values:
    data = []
    for mu in means:
        # To achieve varying mean with fixed p, vary n_trials per group
        trials = int(mu / p)
        samples = np.random.binomial(trials, p, size=group_size)
        data.append((samples.mean(), samples.var()))
    means_, vars_ = zip(*data)
    # Add epsilon to ensure positive values
    log_means = np.log10(np.maximum(means_, EPSILON))
    log_vars = np.log10(np.maximum(vars_, EPSILON))
    slope, intercept, r_value, _, _ = linregress(log_means, log_vars)
    b_values_binom.append(slope)
    r2_values_binom.append(r_value**2)

for r in dispersion_values:
    data = []
    for mu in means:
        p = r / (mu + r)
        samples = nbinom.rvs(r, p, size=group_size)
        data.append((samples.mean(), samples.var()))
    means_, vars_ = zip(*data)
    # Add epsilon to ensure positive values
    log_means = np.log10(np.maximum(means_, EPSILON))
    log_vars = np.log10(np.maximum(vars_, EPSILON))
    slope, intercept, r_value, _, _ = linregress(log_means, log_vars)
    b_values_nb.append(slope)
    r2_values_nb.append(r_value**2)

for sigma in sigma_values:
    data = []
    for mu in means:
        scale = np.exp(np.log(mu) - 0.5 * sigma**2)
        samples = lognorm.rvs(sigma, scale=scale, size=group_size)
        data.append((samples.mean(), samples.var()))
    means_, vars_ = zip(*data)
    # Add epsilon to ensure positive values
    log_means = np.log10(np.maximum(means_, EPSILON))
    log_vars = np.log10(np.maximum(vars_, EPSILON))
    slope, intercept, r_value, _, _ = linregress(log_means, log_vars)
    b_values_ln.append(slope)
    r2_values_ln.append(r_value**2)

for std in std_devs:
    data = []
    for mu in means:
        samples = np.random.normal(loc=mu, scale=std, size=group_size)
        data.append((samples.mean(), samples.var()))
    means_, vars_ = zip(*data)
    # Add epsilon to ensure positive values
    log_means = np.log10(np.maximum(means_, EPSILON))
    log_vars = np.log10(np.maximum(vars_, EPSILON))
    slope, intercept, r_value, _, _ = linregress(log_means, log_vars)
    b_values_gauss.append(slope)
    r2_values_gauss.append(r_value**2)

# Updated plot with Gaussian included
plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.scatter(b_values_nb, r2_values_nb, color="green", label="Negative Binomial")
plt.xlabel("Taylor's Law Exponent b")
plt.ylabel("R² of Log-Log Fit")
plt.title("Negative Binomial")
plt.grid(True)

plt.subplot(1, 3, 2)
plt.scatter(b_values_ln, r2_values_ln, color="red", label="Log-Normal")
plt.xlabel("Taylor's Law Exponent b")
plt.ylabel("R² of Log-Log Fit")
plt.title("Log-Normal")
plt.grid(True)

plt.subplot(1, 3, 3)
plt.scatter(b_values_gauss, r2_values_gauss, color="blue", label="Gaussian")
plt.xlabel("Taylor's Law Exponent b")
plt.ylabel("R² of Log-Log Fit")
plt.title("Gaussian")
plt.grid(True)

plt.suptitle("Correlation Between Taylor's Law Exponent (b) and R² Across Distributions")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('distributions_comparison.png', dpi=300)
print("Saved plot to distributions_comparison.png")

plt.figure(figsize=(7, 5))
plt.scatter(b_values_binom, r2_values_binom, color="purple", label="Binomial")
plt.xlabel("Taylor's Law Exponent b")
plt.ylabel("R² of Log-Log Fit")
plt.title("Binomial Distribution (Underdispersion)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('binomial_distribution.png', dpi=300)
print("Saved plot to binomial_distribution.png")
