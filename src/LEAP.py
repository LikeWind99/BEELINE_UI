import pandas as pd
import numpy as np


class LEAP:
    def __init__(self, pseudotimeData, exprmatxData):
        self.pseudotimeData = pseudotimeData
        self.exprmatxData = exprmatxData

    def MAC_lags(self,
                 pseudotimeData,
                 exprmatxData,
                 max_lag_prop=1 / 3,
                 symmetric=True):
        ###############################################
        # 数据读入
        ###############################################
        # pseudotimeData = pd.read_csv(self.pseudotime_dir, header=0, names=[
        #                           'cell', 'PseudoTime1', 'PseudoTime2'])
        # pseudotimeData = pseudotimeData.fillna(0)
        # exprmatxData = pd.read_csv('ExpressionData.csv', index_col=0, header=0)

        pseudotimeData['PseudoTime'] = pseudotimeData['PseudoTime1'] + \
            pseudotimeData['PseudoTime2']

        ###############################################
        # 数据预处理
        ###############################################
        pseudotimeData.sort_values(by=['PseudoTime', 'cell'],
                                   ascending=True,
                                   inplace=True)

        # extract pesudotime data
        cells = list(pseudotimeData['cell'])
        # sort cell by pseudotime in the data
        expr_Matx = exprmatxData[cells]
        gene_names = list(expr_Matx.index)

        ###############################################
        # 初始化数据
        ###############################################
        num_genes, num_time_int = expr_Matx.shape
        max_lag = np.floor(num_time_int * max_lag_prop)
        window = num_time_int - max_lag

        ###############################################
        # 标准化数据 -> 计算相关系数
        ###############################################
        # data_0 : shape(num_genes, window)
        data_0 = expr_Matx.iloc[:, 0:window]
        cent_0 = (data_0 - data_0.mean(axis=1)).T
        # cor_0 : shape(num_genes, num_genes)
        cor_0 = cent_0.corr(method='pearson')

        # rowsumx_0, rowsumx2_0 : shape(num_genes, 1)
        rowsumx_0 = np.sum(cor_0, axis=0)
        rowsumx2_0 = np.sum(np.power(cor_0, 2), axis=0)

        ###############################################
        # 计算滞后的相关系数
        ###############################################
        max_lag_cor = cor_0
        lag_max = np.zeros((num_genes, num_genes))

        # lag_max_record = np.zeros([num_genes, num_genes])

        for i in range(max_lag):
            ###############################################
            # 获得矩阵并归一化
            ###############################################
            data_lag = expr_Matx.iloc[:, (i + 1):(window + i + 1)]

            # cent_lag : shape(num_genes, window)
            cent_lag = (data_lag - data_lag.mean(axis=1)).T

            # rowsumx_lag, rowsumx2_lag : shape(num_genes, 1)
            rowsumx_lag = np.sum(cent_lag, axis=1)
            rowsumx2_lag = np.sum(np.power(cent_lag, 2), axis=1)

            cor_lag = (cent_0.T.dot(cent_lag) -
                       1 / window * rowsumx_0.T.dot(rowsumx_lag)) / (np.sqrt(
                           rowsumx2_0 - np.power(rowsumx_0, 2) / window).T.dot(
                               np.sqrt(rowsumx2_lag -
                                       np.power(rowsumx_lag, 2) / window)))

            # Keep cor only if greater than previous lag cor
            max_lag_cor[max_lag_cor.abs() < cor_lag.abs()] = cor_lag[
                max_lag_cor.abs() < cor_lag.abs()]
            lag_max[max_lag_cor.abs() < cor_lag.abs()] = i

        MACs_greatest = max_lag_cor
        lag_greatest = lag_max

        if symmetric:
            # max_lag_cor : shape(num_genes, num_genes)
            max_lag_cor_t = max_lag_cor.T
            lag_max_t = lag_max.T
            # 获得横坐标
            coor_x, coor_y = np.where(max_lag_cor.abs() >= max_lag_cor_t.abs())

            MACs_greatest.iloc[coor_x, coor_y] = max_lag_cor.iloc[coor_x,
                                                                  coor_y]
            lag_greatest.iloc[coor_x, coor_y] = lag_max.iloc[coor_x, coor_y]
            MACs_greatest.iloc[-coor_x, -coor_y] = max_lag_cor_t.iloc[-coor_x,
                                                                      -coor_y]
            lag_greatest.iloc[-coor_x, -coor_y] = lag_max_t.iloc[-coor_x,
                                                                 -coor_y]
        np.fill_diagonal(MACs_greatest, np.nan)
        np.fill_diagonal(lag_greatest, np.nan)
        result = pd.DataFrame(MACs_greatest,
                              index=gene_names,
                              columns=gene_names)
        return result


if __name__ == '__main__':
    workflow = LEAP('PseudoTime.csv', 'ExpressionData.csv')
    interaction_matx = workflow.MAC_lags()
