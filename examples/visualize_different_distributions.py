import numpy as np
import matplotlib.pyplot as plt
import math

# Assume the normal (mean) operation time is x
x = 10  # for example, 10 time units

num_samples = 10000

# 1. Normal Distribution
mu_normal = x
sigma_normal = 0.1 * x  # assume 10% uncertainty
samples_normal = np.random.normal(mu_normal, sigma_normal, num_samples)
samples_normal = samples_normal[samples_normal >= 0]  # remove negatives if necessary

# 2. Lognormal Distribution
sigma_lognormal = 0.3  # chosen uncertainty on log scale
mu_lognormal = np.log(x) - (sigma_lognormal**2) / 2
samples_lognormal = np.random.lognormal(mean=mu_lognormal, sigma=sigma_lognormal, size=num_samples)

# 3. Exponential Distribution
scale_exponential = x  # mean equals the scale parameter
samples_exponential = np.random.exponential(scale=scale_exponential, size=num_samples)

# 4. Gamma Distribution
k_gamma = 2  # chosen shape parameter
theta_gamma = x / k_gamma  # scale parameter computed to give mean = x
samples_gamma = np.random.gamma(shape=k_gamma, scale=theta_gamma, size=num_samples)

# 5. Weibull Distribution
k_weibull = 1.5  # chosen shape parameter
lambda_weibull = x / math.gamma(1 + 1/k_weibull)  # scale parameter computed from mean
samples_weibull = np.random.weibull(a=k_weibull, size=num_samples) * lambda_weibull

# Plot the histograms for each distribution
plt.figure(figsize=(12, 8))
plt.hist(samples_normal, bins=50, alpha=0.5, label='Normal (truncated)')
plt.hist(samples_lognormal, bins=50, alpha=0.5, label='Lognormal')
plt.hist(samples_exponential, bins=50, alpha=0.5, label='Exponential')
plt.hist(samples_gamma, bins=50, alpha=0.5, label='Gamma')
plt.hist(samples_weibull, bins=50, alpha=0.5, label='Weibull')
plt.xlabel('Operation Time')
plt.ylabel('Frequency')
plt.title('Operation Time Distributions with Mean = x')
plt.legend()
plt.show()
