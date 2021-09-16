import pandas as pd
import numpy as np
import numpy
from sklearn import linear_model


class SCODE:
    def __init__(self, pseudotimeData, exprmatxData, pnum=4, epochs=100):
        self.pseudotimeData = pseudotimeData
        self.exprmatxData = exprmatxData
        self.tfnum, self.cnum = exprmatxData.shape
        self.pnum = pnum
        self.epochs = epochs

    def sampleZ(self, outIter, inIter, x, y, z):
        for i in range(outIter):
            for j in range(inIter):
                x[i,
                  j] = np.exp(y[i] * z[j]) + np.random.uniform(-0.001, 0.001)

    # SCODE核心函数
    def SCODE_Run(self, exprmatxData, pseudotimeData):

        ###############################################
        # 数据初始化
        ###############################################
        exprmatxData = exprmatxData.values
        W = np.zeros((self.tfnum, self.pnum), dtype=np.float32)
        Z = np.zeros((self.pnum, self.cnum), dtype=np.float32)
        WZ = np.empty((self.tfnum, self.cnum), dtype=np.float32)

        maxB = 2.0
        minB = -10.0
        new_B = np.zeros(self.tfnum, dtype=np.float32)
        old_B = np.zeros(self.tfnum, dtype=np.float32)
        RSS = np.float('inf')

        pseudotimeData['PseudoTime'] = pseudotimeData['PseudoTime1'] + \
            pseudotimeData['PseudoTime2']
        pstime = pseudotimeData['PseudoTime'] / np.max(
            pseudotimeData['PseudoTime'])

        ###############################################
        # 参数初始化
        ###############################################
        for i in range(1, self.pnum):
            old_B[i] = new_B[i] = np.random.uniform(minB, maxB)

        ###############################################
        # 迭代拟合
        ###############################################
        for epoch in range(self.epochs):
            target = np.floor(np.random.uniform(1, self.pnum + 1))
            new_B[target] = np.random.uniform(minB, maxB)
            if (epoch == self.epochs):
                new_B = old_B
            self.sampleZ(self.pnum, self.cnum, Z, new_B, pstime)
            for i in range(self.tfnum):
                reg = linear_model.LinearRegression()
                reg.fit(numpy.transpose(Z), exprmatxData[i, ])
                for j in range(self.pnum):
                    W[i, j] = reg.coef_[j]
                WZ[i, :] = W[i, :] @ Z

            tmp_RSS = np.sum(np.power((exprmatxData - WZ), 2))
            if tmp_RSS < RSS:
                RSS = tmp_RSS
                old_B[target] = new_B[target]
            else:
                new_B[target] = old_B[target]

            B = np.zeros((self.pnum, self.pnum), dtype=np.float32)
            np.fill_diagonal(B, new_B)

            invW = np.linalg.pinv(W)
            A = W @ B @ invW
            return A


if __name__ == '__main__':
    workflow = SCODE('E:/PseudoTime.csv', 'E:/ExpressionData.csv', 19, 4, 2000,
                     100)
    interaction_matx = workflow.SCODE_Run()
