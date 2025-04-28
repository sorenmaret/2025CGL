import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, nbinom, lognorm
from sklearn.linear_model import LinearRegression

# Settings
np.random.seed(42)
n_groups = 50
group_size = 1000

# Create different means for groups
means = np.linspace(1, 100, n_groups)

# Prepare containers
taylor_data = {
    "Poisson": [],
    "Negative Binomial": [],
    "Log-Normal": []
}

# Generate data for each distribution type
for mu in means:
    # Poisson: Var = mu
    samples_poisson = poisson.rvs(mu, size=group_size)
    taylor_data["Poisson"].append((samples_poisson.mean(), samples_poisson.var()))

    # Negative Binomial: Var > mu, parameterized by mean and dispersion
    r = 10  # dispersion parameter (smaller r = more overdispersion)
    p = r / (mu + r)
    samples_nb = nbinom.rvs(r, p, size=group_size)
    taylor_data["Negative Binomial"].append((samples_nb.mean(), samples_nb.var()))

    # Log-normal: highly skewed, multiplicative processes
    sigma = 1.0  # spread/shape of the lognormal distribution
    scale = np.exp(np.log(mu) - 0.5 * sigma**2)  # so that E[X] ≈ mu
    samples_lognorm = lognorm.rvs(sigma, scale=scale, size=group_size)
    taylor_data["Log-Normal"].append((samples_lognorm.mean(), samples_lognorm.var()))

# Function to fit Taylor's Law and get exponent b
def fit_taylor_law(data):
    means, variances = zip(*data)
    X = np.log10(means).reshape(-1, 1)
    y = np.log10(variances)
    reg = LinearRegression().fit(X, y)
    return reg.coef_[0], means, variances  # return b, and data for plotting

# Fit and plot
plt.figure(figsize=(10, 6))
colors = {"Poisson": "blue", "Negative Binomial": "green", "Log-Normal": "red"}

for dist, data in taylor_data.items():
    b, means, vars_ = fit_taylor_law(data)
    plt.scatter(means, vars_, label=f"{dist} (b = {b:.2f})", alpha=0.7, color=colors[dist])

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Mean")
plt.ylabel("Variance")
plt.title("Taylor's Law Applied to Different Distributions")
plt.legend()
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.show()

