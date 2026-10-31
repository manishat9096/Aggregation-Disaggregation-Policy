import numpy as np
import scipy.stats as stats

class CapacityMarginal:
    def __init__(self, mu, sd, asset_type='ordinary', p_outage=0.05):
        if sd <= 0 or asset_type not in ('ordinary', 'outage') or not 0 <= p_outage < 1:
            raise ValueError('Invalid marginal parameters')
        self.mu, self.sd, self.asset_type = float(mu), float(sd), asset_type
        self.p_outage = float(p_outage)
        self.positive = stats.truncnorm(-self.mu/self.sd, np.inf, loc=self.mu, scale=self.sd)
        self.prob_zero = (float(stats.norm.cdf(-self.mu/self.sd))
                          if asset_type == 'ordinary' else self.p_outage)

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        if self.asset_type == 'ordinary':
            value = stats.norm.cdf(x, loc=self.mu, scale=self.sd)
        else:
            value = self.p_outage + (1-self.p_outage)*self.positive.cdf(x)
        return np.where(x < 0, 0., value)

    def ppf(self, u):
        u = np.asarray(u, dtype=float)
        if np.any((u < 0) | (u > 1) | ~np.isfinite(u)):
            raise ValueError('Quantile probabilities must be in [0, 1]')
        if self.asset_type == 'ordinary':
            return np.maximum(0., stats.norm.ppf(u, loc=self.mu, scale=self.sd))
        v = np.clip((u-self.p_outage)/(1-self.p_outage), 0., 1.)
        return np.where(u <= self.p_outage, 0., np.maximum(0., self.positive.ppf(v)))

    def sample(self, n, rng):
        return self.ppf(rng.random(n))

    def failure(self, c):
        # P(X < c) = F(c-); the only atom is at zero.
        c = np.asarray(c, dtype=float)
        return np.where(c <= 0, 0., self.cdf(c))