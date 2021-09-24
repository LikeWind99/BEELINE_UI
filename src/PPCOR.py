import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
from scipy.stats import norm


class PCOR:
    def __init__(self, x, method='pearson'):
        self.method = method
        self.x = x

    def Pcor(self, x, method):
        if not isinstance(x, pd.DataFrame):
            x = pd.DataFrame(x)
        ###############################################
        # 数据初始化
        ###############################################
        n, gp = x.shape[0], x.shape[1] - 2
        cvx = x.cov()

        icvx = pd.DataFrame(np.linalg.inv(cvx))
        pcor = -icvx.corr(method=method)
        np.fill_diagonal(pcor.values, 1)

        if method == 'kendall':
            statistic = pcor / \
                np.sqrt(2 * 2 * (n - gp) + 5) / (9 * (n - gp) * (n - 1 - gp))
            p_value = 2 * norm.cdf(-np.abs(statistic))
        else:
            statistic = pcor * np.sqrt((n - 2 - gp) / (1 - np.power(pcor, 2)))
            p_value = 2 * norm.cdf(-np.abs(statistic))
            # print(p_value)

        np.fill_diagonal(statistic.values, 0)
        np.fill_diagonal(p_value, 0)
        # df = {'estimate': pcor,
        #       'p_value': p_value,
        #       'statistic': statistic,
        #       'n': n,
        #       'gp': gp,
        #       'method': method}
        return p_value

    def spcor(self, x, method):
        n, gp = x.shape[0], x.shape[1] - 2
        x = pd.DataFrame(x)
        cvx = x.cov()

        icvx = pd.DataFrame(np.linalg.inv(cvx))
        spcor = -icvx.corr(method=method)/np.sqrt(np.diag(cvx)) / \
            np.sqrt(np.abs(np.diag(icvx) - ((icvx ** 2).T / np.diag(icvx)).T))
        np.fill_diagonal(spcor.values, 1)

        if method == 'kendall':
            statistic = spcor / \
                np.sqrt(2 * (2 * (n - gp) + 5) / (9 * (n - gp) * (n - 1 - gp)))
            p_value = 2 * norm.cdf(-np.abs(statistic))
        else:
            statistic = spcor * \
                np.sqrt((n - 2 - gp) / (1 - np.power(spcor, 2)))
            p_value = 2 * norm.cdf(-np.abs(statistic))
        np.fill_diagonal(statistic.values, 0)
        np.fill_diagonal(p_value, 0)

        # df = {'estimate': spcor,
        #       'p_value': p_value,
        #       'statistic': statistic,
        #       'n': n,
        #       'gp': gp,
        #       'method': method}
        # print(x.columns.values)
        df = pd.DataFrame(p_value,columns=x.columns.values, index=x.columns.values)
        # print(df)
        return df


if __name__ == '__main__':
    X = np.random.randn(100, 20)
    sample = PCOR(X)
    print(sample.spcor(X, sample.method).shape)
