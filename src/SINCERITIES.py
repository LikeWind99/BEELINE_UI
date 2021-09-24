from numpy.core.einsumfunc import _parse_possible_contraction
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import RidgeCV


class DATA:
    def __init__(self, filePath):
        self.filePath = filePath

    def upload(self, filePath):
        # 获取原始数据
        wholeData = pd.read_csv(filePath)

        # 获得基因名称
        genes = wholeData.columns.values[:-1]

        # 获得基因数量
        numGENES = len(genes)

        wholeData.dropna(axis=1, how='any', inplace=True)
        timeLine = wholeData.iloc[:, -1].values

        time = sorted(list(set(timeLine)))
        num_time_points = len(time)

        sortTOTdata = []
        sortTIMELINE = []

        for k in range(num_time_points):
            sortTOTdata.append(wholeData[wholeData.iloc[:, -1] == time[k]])
            sortTIMELINE.append(
                wholeData[wholeData.iloc[:, -1] == time[k]].iloc[:, -1])

        return time, num_time_points, sortTOTdata, sortTIMELINE, numGENES, genes


class sincerities:
    def __init__(self, filePath):
        # self.data = pd.read_csv(filePath)
        self.time, self.num_time_points, self.sortTIMELINE, self.genes, self.numGENES, self.singleCELLdata = DATA(
            filePath).upload(filePath)
        # self.time = time  # list
        # self.num_time_points = num_time_points  # length of time 8
        # self.sortTIMELINE = sortTIMELINE  # df shape:(960, 1)
        # self.genes = genes  # list
        # self.numGENES = numGENES  # length of genes: 45
        # self.singleCELLdata = sortTOTdata  # list 8 each:df shape:(120, 46)

    def sinceri(self, distance=1, SIGN=1):
        single_cell_data = self.singleCELLdata
        time = self.time
        numGENES = self.numGENES
        num_time_points = self.num_time_points

        DISTANCE_matrix = np.zeros(
            (num_time_points - 1, numGENES))  # shape(7, 45)
        totalDATA = single_cell_data[0]
        for ti in range(num_time_points - 1):
            totalDATA = np.append(totalDATA, single_cell_data[ti + 1])
            data_ti = single_cell_data[ti].T

            data_ti_plus1 = single_cell_data[ti + 1].T

            for gi in range(numGENES):
                p1 = data_ti.iloc[gi, :].values
                p2 = data_ti_plus1.iloc[gi, :].values
                if (distance == 1):
                    p1 = p1.flatten()
                    p2 = p2.flatten()
                    test_stat = stats.ks_2samp(p1, p2)
                    DISTANCE_matrix[ti, gi] = test_stat.statistic

        time = np.array(time)
        persudotime = (time[1:] - time[:-1]).reshape(-1, 1)

        deltaT = np.zeros((num_time_points - 1, numGENES))
        deltaT[:, :] = persudotime

        DISTANCE_matrix = DISTANCE_matrix / deltaT

        X_matrix = DISTANCE_matrix[:(num_time_points - 2), :]

        pred_lambda_min = np.zeros((numGENES, numGENES))

        alphas = np.logspace(-10, -0.01, 100)
        lambda_res = []

        for gi in range(numGENES):
            # beta = np.zeros((X_matrix.shape[1], 1))
            Y_vector = DISTANCE_matrix[1:(num_time_points), gi]
            ridge = RidgeCV(alphas=alphas).fit(X_matrix, Y_vector)
            pred_lambda_min[:, gi] = ridge.coef_

            lambda_res.append(ridge.alpha_)

        max_pred = np.amax(pred_lambda_min)
        # print(max_pred)
        adj_matrix = pred_lambda_min / max_pred
        # interactions = np.array(adj_matrix.flatten())
        # print(len(interactions))
        # edges = np.array(["no regulation"] * len(interactions))
        # SIGN = 1
        # if SIGN == 1:
        #     gt_loc = np.where(interactions > 0)
        #     for x in gt_loc:
        #         edges[x] = "activation"
        #     lt_loc = np.where(interactions < 0)
        #     for x in lt_loc:
        #         edges[x] = "repression"

        # else:
        #     ne0_loc = np.where(interactions != 0)
        #     for x in ne0_loc:
        #         edges[x] = "activation/repression"

        # interactions = abs(interactions)
        # targetGENES = list(self.genes) * 45
        # sourceGENES = [val for val in self.genes for i in range(45)]

        # df = pd.DataFrame(columns=('targetGENES', 'sourceGENES', 'Interaction',
        #                            'Edges'))
        # df = df.append(
        #     pd.DataFrame({
        #         'targetGENES': targetGENES,
        #         'sourceGENES': sourceGENES,
        #         'Interaction': interactions,
        #         'Edges': edges
        #     }))
        # df = df.sort_values(axis=0, by='Interaction', ascending=False)
        # df.index = range(1, 2026)
        return adj_matrix


if __name__ == '__main__':
    filePath = './data/THP1_single_cell_data.csv'
    sincerity = sincerities(filePath)
    sincerity.sinceri()
