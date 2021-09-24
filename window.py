import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns

from PyQt5 import uic
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *

import matplotlib
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from src import help
from src import LEAP
from src import PPCOR
from src import SCODE
from src import design
from src import picture
from src import GRNBoost2
from src import SINCERITIES

from scipy.sparse import issparse

import scanpy as sc
from scanpy.plotting._tools.scatterplots import tsne
from scanpy.preprocessing._normalization import normalize_total

from pyscenic.aucell import aucell
from pyscenic.prune import prune2df, df2regulons
from pyscenic.utils import modules_from_adjacencies
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase


matplotlib.use('Qt5Agg')  # 使用Qt后台，以便嵌入Qt
font = FontProperties(
    fname=r"c:\windows\fonts\simsun.ttc", size=10)  # 设置字体用于显示中文


class mainWindow(QtWidgets.QMainWindow, design.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        # 不同版面图片存储的前缀
        self.picName = ['Weighted_graph', 'AUCell', 't-SNE']

        # 语言切换
        self.trans = QTranslator()  # 创建翻译器
        self.ui_plotArea_1 = QStackedLayout(
            self.plotArea_1)  # 第一个版面作图区域创建stack布局
        self.ui_plotArea_2 = QStackedLayout(
            self.plotArea_2)  # 第二个版面作图区域创建stack布局
        self.ui_plotArea_3 = QStackedLayout(
            self.plotArea_3)  # 第三个版面作图区域创建stack布局

        # 左侧版面切换按钮绑定
        self.btn1.clicked.connect(self.btn1_fun)  # 切换为第一个版面
        self.btn2.clicked.connect(self.btn2_fun)  # 切换为第二个版面
        self.btn3.clicked.connect(self.btn3_fun)  # 切换为第三个版面

        # 上方图标函数绑定
        self.actionChinese.triggered.connect(
            self.English2Chinese)  # 英文转换为中文(实际上为任意语言转换为中文)
        self.actionEnglish.triggered.connect(
            self.Chinese2English)  # 中文转换为英文(实际上为任意语言转换为英文)
        self.actionexpression_data.triggered.connect(
            self.loadExpressData)  # 表达数据导入
        self.actionpseudotime_data.triggered.connect(
            self.loadPseudoData)  # 伪时间数据导入
        self.actiondocument.triggered.connect(self.loadDoc)  # 加载帮助文档
        self.actioncsv.triggered.connect(self.exportCSV)  # 关系矩阵保存为csv
        self.actionexcel.triggered.connect(self.exportExcel)  # 关系矩阵保存为excel格式
        self.actionPNG.triggered.connect(lambda: self.exportPicPNG(
            self.picName[self.stackedWidget.currentIndex()]))  # 图片保存为png格式
        self.actionJPG.triggered.connect(lambda: self.exportPicJPG(
            self.picName[self.stackedWidget.currentIndex()]))  # 图片保存为jpg格式
        self.actionSVG.triggered.connect(lambda: self.exportPicSVG(
            self.picName[self.stackedWidget.currentIndex()]))  # 图片保存为svg格式

        # 版面一按钮绑定
        self.rbtn1_2.clicked.connect(self.click_1)  # run按钮绑定算法运行线程
        self.rbtn2_2.clicked.connect(
            self.network_visualization)  # plot按钮绑定作图事件

        # 版面二按钮绑定
        self.rbtn1_3.clicked.connect(self.click_2)
        self.rbtn2_3.clicked.connect(self.aucellHeatmap)

        # 版面三按钮绑定
        self.btn4_3.clicked.connect(self.click_3)
        self.btn4_2.clicked.connect(self.tsne_plot)

        # 设置标准输出和标准错误输出的位置为界面的信息框
        sys.stdout = EmittingStream(textWritten=self.outputWritten)
        sys.stderr = EmittingStream(textWritten=self.outputWritten)

    # 在框中输出算法运行过程
    def outputWritten(self, text):
        cursor = self.outputArea.textCursor()  # 获取cursor
        cursor.movePosition(QTextCursor.End)  # 将cursor移动至当前文字的末尾
        cursor.insertHtml(text)  # 在文字后插入输出信息
        self.outputArea.setTextCursor(cursor)  # 设置cursor(必须设置)
        self.outputArea.ensureCursorVisible()  # 确保信息可见

    # 算法运行事件函数
    def click_1(self):
        self.thread_1 = calcRelationMat(self)  # 创建实例
        self.thread_1.start()  # 子线程运行

    # 创建RcisTarget子线程并运行
    def click_2(self):
        self.thread_2 = RcisTarget(self)
        self.thread_2.start()

    # 创建t-SNE子线程并运行
    def click_3(self):
        self.btn4_2.setEnabled(False)  # 禁用作图按钮
        self.thread_3 = t_sne(self)
        self.thread_3._signal.connect(self.set_btn)  # 恢复作图按钮
        self.thread_3.start()

    # 设置按钮状态
    def set_btn(self):
        self.btn4_2.setEnabled(True)

    # 切换为第一版面
    def btn1_fun(self):
        self.stackedWidget.setCurrentIndex(0)  # 将当前版面设为第一版面

    # 切换为第二版面
    def btn2_fun(self):
        self.stackedWidget.setCurrentIndex(1)  # 将当前版面设为第二版面

    # 切换为第三版面
    def btn3_fun(self):
        self.stackedWidget.setCurrentIndex(2)  # 将当前版面设为第三版面

    # 中文转换成英文
    def Chinese2English(self):
        self.trans.load('en')  # 加载英文文件
        _app = QApplication.instance()  # 获取实例
        _app.installTranslator(self.trans)  # 翻译
        self.retranslateUi(self)  # 界面重构

    # 英语转换成中文
    def English2Chinese(self):
        self.trans.load('zh_cn')  # 加载简体中文文件
        _app = QApplication.instance()  # 获取实例
        _app.installTranslator(self.trans)  # 翻译
        self.retranslateUi(self)  # 界面重构

    # 推断网络可视化
    def network_visualization(self, threshold=0.75, remapping=False):
        '''
            @param: threshold    阈值，据此判断可能存在的关系
            @param: remapping    重排序(实际上并没有使用(一.一||)，可删除)
        '''
        # 直接判断
        global grnfig
        if not 'relationMatrix' in globals():
            QMessageBox.about(
                self,
                'Warning',
                'please run the algorithm or wait patiently for the algorithm to finish',
            )
            return
        inter_matx = relationMatrix  # 创建局部变量，避免对全局变量修改，以防后面使用
        util = utils()
        inter_matx = util.isDataFrame(inter_matx)  # 确保输入为DataFrame

        GRN = nx.Graph()  # 创建新图
        grnfig, ax = plt.subplots(1, 1)  # 创建画布
        genes = list(inter_matx.columns)  # 获取基因名
        # 归一化
        inter_matx = (inter_matx - np.min(inter_matx)) / (np.max(inter_matx) -
                                                          np.min(inter_matx))

        percenile_val = np.percentile(inter_matx, 100 * threshold)  # 计算百分比
        # 设置cmap
        cmap = colors.ListedColormap([
            '#FF300E', '#FF6B01', '#FFE500', '#7FF200', '#40F200', '#04CD65',
            '#00A3A2', '#006AE0', '#4D25FF', '#7A14F3'
        ])
        # 设置cNorm
        cNorm = colors.BoundaryNorm(
            [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], cmap.N)
        # 从给定的colormap返回RGBA颜色之前使用数据归一化化
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
        # 权重列表
        weightList = []
        # 判断需要画出的边和点
        for i in range(inter_matx.shape[0]):
            for j in range(i + 1, inter_matx.shape[1]):
                if inter_matx.iloc[i, j] > percenile_val:
                    GRN.add_edge(genes[i], genes[j])
                    weightList.append(inter_matx.iloc[i, j])
        if remapping:
            weightList = (weightList - np.min(weightList)) / (
                np.max(weightList) - np.min(weightList))
            colorList = list(
                map(lambda weight: scalarMap.to_rgba(weight), weightList))
        else:
            colorList = list(
                map(lambda weight: scalarMap.to_rgba(weight), weightList))
        # 创建布局
        pos = nx.spring_layout(GRN)
        # 画点
        nx.draw_networkx_nodes(GRN,
                               pos,
                               ax=ax,
                               node_size=900,
                               node_color='#99ccff',
                               alpha=0.7)
        # 画边
        nx.draw_networkx_edges(GRN,
                               pos,
                               ax=ax,
                               edgelist=GRN.edges,
                               width=0.8,
                               edge_color=colorList)
        # 画注释
        nx.draw_networkx_labels(GRN,
                                pos,
                                ax=ax,
                                font_size=6,
                                font_family='sans-serif',
                                font_weight='heavy')
        plt.colorbar(scalarMap)  # 调出比色板
        plt.axis('off')  # 关闭坐标轴
        fig = FigureCanvasQTAgg(grnfig)  # 转化为Qt可用格式
        try:
            self.ui_plotArea_1.removeWidget(self.ui_plotArea_1.currentWidget())
            self.ui_plotArea_1.addWidget(fig)  # 加入到界面中
        except:
            self.ui_plotArea_1.addWidget(fig)  # 加入到界面中
        QApplication.processEvents()

    # 加载表达数据
    def loadExpressData(self):
        # 将expression_data和expressionFilePath定义为全局变量
        global expression_data, expressionFilePath

        # 使用try-except捕捉异常并反馈
        try:
            expressionFilePath = QFileDialog.getOpenFileName(
                self, '选择文件', '', '*.csv ')
            # print(expressionFilePath)  # 检查路径是否正确

            util = utils()
            expression_data = util.loadData(expressionFilePath[0]).T
            # print(expression_data)  # 检查数据类型是否为DataFrame
        except OSError:
            print("Please <strong>select the correct file path</strong>!<br>")  # 打印异常

    # 加载伪时间数据
    def loadPseudoData(self):
        # 将pseudotime_data和pseudotimeFilePath定义为全局变量
        global pseudotime_data, pseudotimeFilePath

        # 使用try-except捕捉异常并反馈
        try:
            pseudotimeFilePath = QFileDialog.getOpenFileName(
                self, '选择文件', '', '*.csv ')
            # print(filePath)
            util = utils()
            pseudotime_data = util.loadData(pseudotimeFilePath[0])

        except OSError:
            print("Please <strong>select the correct file path</strong>!<br>")  # 打印异常

    # 实现重名存储
    def findFileName(self, filename, fileType):
        '''
            查找当前目录下重名文件的数量
            @param: filename    重名文件的名称
            @param: fileType    重名文件的类型(考虑到后面需要多次调用，因此写入函数封装)
            @return: count      重名文件的数量
        '''
        match = r'%s.*\.%s' % (filename, fileType)
        result = [
            re.findall(match, file)[0] for file in os.listdir('./figures')
            if len(re.findall(match, file))
        ]
        count = len(result)

        return count

    # 图片导出为SVG格式
    def exportPicSVG(self, filename):
        # print('I am clicked!')  # 测试按钮是否绑定成功
        # 使用try-except捕捉异常
        try:

            count = self.findFileName(filename, 'svg')  # 查找目前存在多少张重名图片
            idx = self.stackedWidget.currentIndex()
            if idx == 0:
                if count:
                    grnfig.savefig(
                        f"./figures/{filename}_{count}.svg", dpi=300)  # 图片重名存储
                else:
                    grnfig.savefig(f"./figures/{filename}.svg", dpi=300)
            elif idx == 1:
                if count:
                    aucellfig.savefig(
                        f"./figures/{filename}_{count}.svg", dpi=300)  # 图片重名存储
                else:
                    aucellfig.savefig(f"./figures/{filename}.svg", dpi=300)
            else:
                pass

        except Exception as e:
            print(e)  # 捕捉异常并打印，实际上永远不会执行(除非出现不可预知的错误~)

    # 图片导出为PNG格式
    def exportPicPNG(self, filename):
        # print('I am clicked!')  # 测试按钮是否绑定成功
        # 使用try-except捕捉异常
        try:
            count = self.findFileName(filename, 'png')  # 查找目前存在多少张重名图片
            idx = self.stackedWidget.currentIndex()
            if idx == 0:
                if count:
                    grnfig.savefig(
                        f"./figures/{filename}_{count}.png", dpi=300)  # 图片重名存储
                else:
                    grnfig.savefig(f"./figures/{filename}.png", dpi=300)
            elif idx == 1:
                if count:
                    aucellfig.savefig(
                        f"./figures/{filename}_{count}.png", dpi=300)  # 图片重名存储
                else:
                    aucellfig.savefig(f"./figures/{filename}.png", dpi=300)
        except Exception as e:
            print(e)  # 捕捉异常并打印，实际上永远不会执行(除非出现不可预知的错误~)

    # 图片导出为JPG格式
    def exportPicJPG(self, filename):
        # print('I am clicked!')  # 测试按钮是否绑定成功
        # 使用try-except捕捉异常
        try:
            count = self.findFileName(filename, 'jpg')  # 查找目前存在多少张重名图片
            idx = self.stackedWidget.currentIndex()
            if idx == 0:
                if count:
                    grnfig.savefig(
                        f"./figures/{filename}_{count}.jpg", dpi=300)  # 图片重名存储
                else:
                    grnfig.savefig(f"./figures/{filename}.jpg", dpi=300)
            elif idx == 1:
                if count:
                    aucellfig.savefig(
                        f"./figures/{filename}_{count}.jpg", dpi=300)  # 图片重名存储
                else:
                    aucellfig.savefig(f"./figures/{filename}.jpg", dpi=300)
        except Exception as e:
            print(e)  # 捕捉异常并打印，实际上永远不会执行(除非出现不可预知的错误~)

    # 关系矩阵导出为CSV格式
    def exportCSV(self):
        try:
            pd.DataFrame(relationMatrix).to_csv('./relationMatrix.csv')
        except NameError:
            print(
                "Please <strong>run the algorithm or wait patiently</strong> for the algorithm to complete!<br>")

    # 关系矩阵导出为Excel格式
    def exportExcel(self):
        try:
            pd.DataFrame(relationMatrix).to_excel('./relationMatrix.xls')
        except NameError:
            print(
                "Please <strong>run the algorithm or wait patiently</strong> for the algorithm to complete!<br>")

    # 创建新窗口加载帮助文档
    def loadDoc(self):
        self.child = childWindow()
        self.child.show()

    # 绘制AUC热图
    def aucellHeatmap(self):
        # 使用try-except结构处理异常
        global aucellfig
        try:
            if 'aucMtx' in globals():
                # 调用seaborn中的聚类函数作图，并保存画布
                g = sns.clustermap(aucMtx, figsize=(8, 8))
                # 将画布转换成为Qt控件
                aucellfig = g.fig
                fig = FigureCanvasQTAgg(g.fig)
                # 将控件加入到指定区域
                try:
                    self.ui_plotArea_2.removeWidget(
                        self.ui_plotArea_2.currentWidget())
                    self.ui_plotArea_2.addWidget(fig)  # 加入到界面中
                except:
                    self.ui_plotArea_1.addWidget(fig)  # 加入到界面中
            else:
                # 调用seaborn中的聚类函数作图，并保存画布
                g = sns.clustermap(expression_data, figsize=(8, 8))
                # 将画布转换成为Qt控件
                aucellfig = g.fig
                fig = FigureCanvasQTAgg(g.fig)
                # 将控件加入到指定区域
                try:
                    self.ui_plotArea_2.removeWidget(
                        self.ui_plotArea_2.currentWidget())
                    self.ui_plotArea_2.addWidget(fig)  # 加入到界面中
                except:
                    self.ui_plotArea_2.addWidget(fig)  # 加入到界面中
        except NameError:
            # 如果报错弹窗警告
            QMessageBox.about(
                self,
                'Error',
                'Please load expression data!',
            )

    # 画t-SNE图
    def tsne_plot(self):
        try:
            lable = QtWidgets.QLabel()
            fig = QPixmap('./figures/tsne.png')
            lable.setScaledContents(True)
            lable.setPixmap(fig)
            try:
                self.ui_plotArea_3.removeWidget(
                    self.ui_plotArea_3.currentWidget())
                self.ui_plotArea_3.addWidget(lable)  # 加入到界面中
            except:
                self.ui_plotArea_3.addWidget(lable)  # 加入到界面中

        except NameError:
            # 如果报错弹窗警告
            QMessageBox.about(
                self,
                'Error',
                'Please load expression data!',
            )
        except:
            print("Something wrong!")


class calcRelationMat(QThread):
    '''
        算法选择并运行
    '''
    qmut = QMutex()

    def __init__(self, mainwindow):
        super(calcRelationMat, self).__init__()
        self.mainwindow = mainwindow
        # 底部算法选择按钮
        self.comboBox = self.mainwindow.comboBox_2

    def run(self):
        '''
            共包含五种算法
            @PPCOR          需要输入expression data
            @GRNBoost2      需要输入expression data
            @LEAP           需要输入expression data和pseudotime data
            @SCODE          需要输入expression data和pseudotime data
            @SINCERITIES    需要输入expression data和pseudotime data
            @return         返回推断出来的关系网络
        '''
        global relationMatrix  # 将relationMatrix设置为全局变量，以便主线程画图使用
        self.qmut.lock()  # 线程加锁
        funcName = self.comboBox.currentText()  # 获得当前选择的算法名称

        if funcName == 'PPCOR':
            # 使用try-except处理异常，如果未加载数据则会打印对应信息
            try:
                sample = PPCOR.PCOR(expression_data)  # 创建实例
                print('Start, please wait a moment...<br>')
                relationMatrix = sample.spcor(sample.x, sample.method)  # 算法运行
                print(
                    'End, you can click plot button to visualize the result!<br>')

            except NameError:
                print("please <strong>load expression data</strong>!<br>")

        elif funcName == 'GRNBoost2':
            # 使用try-except处理异常，如果未加载数据则会打印对应信息
            try:
                sample = GRNBoost2.GRNBoost2(expression_data)  # 创建实例
                print('Start, please wait a moment...<br>')
                relationMatrix = sample.rungrnboost(sample.expressData)  # 算法运行
                print(
                    'End, you can click plot button to visualize the result!<br>')
            except NameError:
                print("please <strong>load expression data</strong>!<br>")

        elif funcName == 'LEAP':
            # 使用try-except处理异常，如果未加载数据则会打印对应信息
            try:
                sample = LEAP.LEAP(pseudotime_data, expression_data)  # 创建实例
                print('Start, please wait a moment...<br>')
                relationMatrix = sample.MAC_lags()  # 算法运行
                print(
                    'End, you can click plot button to visualize the result!<br>')
            except NameError:
                print(
                    "please <strong>load expression data and pseudotime data</strong>!<br>")

        elif funcName == 'SCODE':
            # 使用try-except处理异常，如果未加载数据则会打印对应信息
            try:
                sample = SCODE.SCODE(pseudotime_data, expression_data)  # 创建实例
                print('Start, please wait a moment...<br>')
                relationMatrix = sample.SCODE_Run(
                    expression_data, pseudotime_data)  # 算法运行
                print(
                    'End, you can click plot button to visualize the result!<br>')
            except NameError:
                print(
                    "please <strong>load expression data and pseudotime data</strong>!<br>")

        else:
            # 使用try-except处理异常，如果未加载数据则会打印对应信息
            try:
                sample = SINCERITIES.sincerities(expressionFilePath[0])  # 创建实例
                print('Start, please wait a moment...<br>')
                relationMatrix = sample.sinceri()  # 算法运行
                print(
                    'End, you can click plot button to visualize the result!<br>')

            except NameError:
                print(
                    "please <strong>load expression data and pseudotime data</strong>!<br>")
            except:
                print("Something wrong!<br>")
        self.qmut.unlock()  # 线程解锁


class RcisTarget(QThread):

    qmut = QMutex()  # 创建线程

    def __init__(self, mainwindow):
        super(RcisTarget, self).__init__()
        self.mainwindow = mainwindow
        self.combox = self.mainwindow.comboBox

    def run(self):
        global regulons, aucMtx
        self.qmut.lock()  # 线程解锁
        speice = self.combox.currentText()  # 获取当前的物种
        dir_path = os.path.dirname(os.path.abspath(__file__))  # 获取当前运行的绝对路径
        if speice == 'mouse':
            # 获得注释文件
            motifAnnotationsFile = os.path.join(
                dir_path, r"data\mouse\motifs-v9-nr.mgi-m0.001-o0.0.tbl")
            # 获得数据库文件
            DbsPath = os.path.join(dir_path, r"data\mouse\*.feather")
        else:
            # 获得注释文件
            motifAnnotationsFile = os.path.join(
                dir_path, r"\data\human\motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
            # 获得数据库文件
            DbsPath = os.path.join(dir_path, r"data\human\*.feather")
        util = utils()
        try:
            adjacencies = util.matirx2Triplet(relationMatrix)  # 将关系矩阵转换为邻接形式
            dbFnames = glob.glob(DbsPath)  # 获取全部数据库文件地址
            # print(dbFnames)
            dbs = [RankingDatabase(fname=fname, name=fname)
                   for fname in dbFnames]  # 创建数据库
            # print("dbs:", dbs)
            modules = list(modules_from_adjacencies(
                adjacencies, expression_data, rho_mask_dropouts=True))
            ###########################################################
            # 出问题的地方，无法获取modules，问题是数据量太小，需要大数据
            ###########################################################
            # print("modules:", modules)
            dfMotifs = prune2df(dbs, modules, motifAnnotationsFile)
            # dfMotifs = pd.DataFrame(dfMotifs)
            regulons = df2regulons(dfMotifs)
            aucMtx = aucell(expression_data, regulons)
        except NameError:
            print(
                'Please load the expression data or wait patiently for the end of the algorithm!<br>')
        except:
            print("Can't find the regulons!<br>")
        self.qmut.unlock()  # 线程解锁


class t_sne(QThread):
    # qmut = QMutex()  # 创建线程
    _signal = pyqtSignal()

    def __init__(self, mainWindow):
        super(t_sne, self).__init__()
        self.mainwindow = mainWindow
        self.combox = self.mainwindow.comboBox_3

    def run(self):
        # self.qmut.lock()  # 线程解锁
        index = int(self.combox.currentText())
        if not 'expression_data' in globals():
            return
        # 获取细胞信息
        cellinfo = pd.DataFrame(
            expression_data.index, index=expression_data.index, columns=['sample_index'])
        # 获取基因信息
        geneinfo = pd.DataFrame(expression_data.columns, index=expression_data.columns,
                                columns=['genes_index'])
        # 数据格式转换
        adata = sc.AnnData(expression_data, obs=cellinfo, var=geneinfo)
        # 标准化
        norm_dict = normalize_total(adata, target_sum=100, inplace=False)
        if issparse(norm_dict['X']):
            mean_percent = norm_dict['X'].mean(axis=0).A1
            top_idx = np.argsort(mean_percent)[::-1][:index]

        else:
            mean_percent = norm_dict['X'].mean(axis=0)
            top_idx = np.argsort(mean_percent)[::-1][:index]
        # 获取表达量最大的
        columns = (
            adata.var_names[top_idx]
        )

        file_path = os.path.join(os.path.abspath(__file__), 'figures/tsne.png')
        if os.path.exists(file_path):
            os.remove(file_path)

        sc.tl.tsne(adata)
        sc.pl.tsne(adata=adata, color=columns.values,
                   show=False, save='.png')

        # self.qmut.unlock()  # 线程解锁
        self._signal.emit()


class EmittingStream(QObject):
    '''
        发送信号
    '''
    textWritten = pyqtSignal(str)  # 定义一个发送str的信号

    # 发射信号
    def write(self, text):
        '''
            @param: text    需要发送的文字信息
        '''
        self.textWritten.emit(str(text))


class childWindow(QtWidgets.QMainWindow, help.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        with open(r'说明文档.txt', 'r', encoding='utf-8') as f:
            text = f.read()  # 读入帮助文档内容
        self.textBrowser.insertHtml(text)  # 插入文档为HTML格式


class utils:
    '''
        数据读入基类
    '''

    def loadData(self, dataPath):
        '''
            文件加载函数
            @param: dataPath      文件的存储路径
            @return: DataFrame    返回一个补全的DataFrame
        '''
        data = pd.read_csv(dataPath, sep=',', index_col=0, header=0)
        data = data.fillna(0)
        return data

    def isDataFrame(self, data):
        '''
            判断输入是否为DataFrame，如果不是则将其转换为DataFrame
            @param: data          传入数据
            @return: DataFrame    返回DataFrame
        '''

        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame(data)

        return data

    def matirx2Triplet(self, matrix):
        '''
            将输入的矩阵转换成为邻接表的形式，最终返回一个DataFrame
            @param: matrix    传入数据，理论上讲可传入任意格式
            @return: df       返回邻接形式的DataFrame
        '''
        df = self.isDataFrame(matrix)
        df = df.to_dict(orient='index')
        df = [[key1, key2, df[key1][key2]]
              for key1 in df.keys() for key2 in df[key1].keys()]

        df = pd.DataFrame(df, columns=['TF', 'target', 'importance'])
        # print(df)
        return df


if __name__ == '__main__':
    # 创建实例
    app = QApplication([])

    # 创建类实例
    stats = mainWindow()

    # 展示ui
    stats.show()

    # 进入事件循环
    app.exec_()
