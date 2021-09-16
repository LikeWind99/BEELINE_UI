from arboreto.algo import grnboost2


class GRNBoost2:
    def __init__(self, expressData):
        self.expressData = expressData

    # ex_matrix：表达矩阵，行为样本，列为基因
    def rungrnboost(self, ex_matrix):

        tf_names = list(ex_matrix.columns)

        # compute the GRN
        network = grnboost2(expression_data=ex_matrix, tf_names=tf_names)
        network.columns = ['Gene1', 'Gene2', 'EdgeWeight']

        # 转换成n*n矩阵，其中行基因调控列基因
        matrix = network.pivot('Gene1', 'Gene2', values='EdgeWeight')
        idx = matrix.columns.union(matrix.index)
        matrix = matrix.reindex(index=idx, columns=idx, fill_value=0)
        matrix = matrix.fillna(0)
        return matrix


if __name__ == '__main__':
    GRNMatrix = GRNBoost2(r'PseudoTime.csv')
    result = GRNMatrix.rungrnboost(GRNMatrix.ex_matrix)
    print(result)
