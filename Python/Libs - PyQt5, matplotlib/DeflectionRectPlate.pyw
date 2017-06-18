import sys, os.path, json, math, configparser
from subprocess import Popen, PIPE, STDOUT
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt, numpy as np
from matplotlib import cm
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QApplication, QMessageBox
from enum import Enum

w_JSON_FILE = 'w.json'
r_JSON_FILE = 'r.json'
INPUT_DATA_INI_FILE = 'input_data.ini'
WINDOW_UI_FILE = 'window.ui'

class PlotType(Enum):
    surface = 1
    contour = 2

def getListFromJSON(fname, key):
    tmp = None
    with open(fname, 'r') as f:
        s = f.read()
        tmp = json.loads(s)[key]
    return tmp

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        window = uic.loadUi(WINDOW_UI_FILE, self)
        self.loadInputData(window, INPUT_DATA_INI_FILE)
        self.deleteJSONFiles(w_JSON_FILE, r_JSON_FILE)
        self.inputData = ()

        window.lineEditNumIters.setDisabled(True)
        window.pushButtonPlot.clicked.connect(self.plot)
        window.pushButtonSaveInputData.clicked.connect(self.saveInputData)
        window.pushButtonDelInputData.clicked \
            .connect(lambda: self.delInputData(window))
        
        window.radioButtonCtConst.clicked \
            .connect(lambda: self.toggleLineEditNumIters('radioButtonCtConst'))
        window.radioButtonCtConstContour \
            .clicked.connect \
            (lambda: self.toggleLineEditNumIters('radioButtonCtConstContour'))
        window.radioButtonCt.clicked \
            .connect(lambda: self.toggleLineEditNumIters('radioButtonCt'))
        window.radioButtonCtContour.clicked \
            .connect(lambda: self.toggleLineEditNumIters('radioButtonCtContour'))
        window.radioButtonKarman.clicked \
            .connect(lambda: self.toggleLineEditNumIters('radioButtonKarman'))
        window.radioButtonKarmanContour.clicked \
            .connect(lambda: self.toggleLineEditNumIters('radioButtonKarmanContour'))
        
        self.show()

    def deleteJSONFiles(self, fname1, fname2):
        if os.path.isfile(fname1):
            os.remove(fname1)
        if os.path.isfile(fname2):
            os.remove(fname2)

    def loadInputData(self, window, fname):
        if os.path.isfile(fname):
            config = configparser.ConfigParser()
            config.read(fname)
            window.lineEdit_n.setText(config.get('Input Data', 'n'))
            window.lineEdit_a.setText(config.get('Input Data', 'a'))
            window.lineEdit_b.setText(config.get('Input Data', 'b'))
            window.lineEdit_h.setText(config.get('Input Data', 'h'))
            window.lineEdit_q0.setText(config.get('Input Data', 'q0'))
            window.lineEdit_E.setText(config.get('Input Data', 'E'))
            window.lineEdit_nu.setText(config.get('Input Data', 'nu'))
            window.lineEditNumIters.setText(config.get('Input Data', 'numIters'))

    def saveInputData(self):
        if self.checkLineEdits():
            n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('')
            config = configparser.ConfigParser()

            config.add_section('Input Data')
            config.set('Input Data', 'n', str(n))
            config.set('Input Data', 'a', str(a))
            config.set('Input Data', 'b', str(b))
            config.set('Input Data', 'h', str(h))
            config.set('Input Data', 'q0', str(q0))
            config.set('Input Data', 'E', str(E))
            config.set('Input Data', 'nu', str(nu))
            config.set('Input Data', 'numIters', str(numIters))

            with open(INPUT_DATA_INI_FILE, 'w') as f:
                config.write(f)

    def delInputData(self, window):
        window.lineEdit_n.setText('')
        window.lineEdit_a.setText('')
        window.lineEdit_b.setText('')
        window.lineEdit_h.setText('')
        window.lineEdit_q0.setText('')
        window.lineEdit_E.setText('')
        window.lineEdit_nu.setText('')

        if os.path.isfile(INPUT_DATA_INI_FILE):
            os.remove(INPUT_DATA_INI_FILE)

    def toggleLineEditNumIters(self, rbName):
        if rbName == 'radioButtonCtConst' \
            or rbName == 'radioButtonCtConstContour':
            window.lineEditNumIters.setDisabled(True)
        else:
            window.lineEditNumIters.setDisabled(False)

    def checkLineEdits(self):
        try:
            int(window.lineEdit_n.text())
            float(window.lineEdit_a.text())
            float(window.lineEdit_b.text())
            float(window.lineEdit_h.text())
            float(window.lineEdit_q0.text())
            float(window.lineEdit_E.text())
            float(window.lineEdit_nu.text())
            int(window.lineEditNumIters.text())
            return True
        except ValueError:
            QMessageBox.warning(self, 'Ошибка','Введены неверные данные')
            return False

    def getValues(self, fname_exe):
        n = int(window.lineEdit_n.text())
        a = float(window.lineEdit_a.text())
        b = float(window.lineEdit_b.text())
        h = float(window.lineEdit_h.text())
        q0 = float(window.lineEdit_q0.text())
        E = float(window.lineEdit_E.text())
        nu = float(window.lineEdit_nu.text())
        numIters = int(window.lineEditNumIters.text())
        return (n, a, b, h, q0, E, nu, numIters, fname_exe)

    def getXY(self, n, a, b):
        hx = a / n
        hy = b / n
        x = np.arange(0, a + hx, hx)
        y = np.arange(0, b + hy, hy)
        x, y = np.meshgrid(x, y)
        return (x, y)

    def plot(self):
        if self.checkLineEdits():
            if window.radioButtonCtConst.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('ct_const.exe')
                self.plotCtConst(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(1), w_JSON_FILE, r_JSON_FILE, fname_exe)
            if window.radioButtonCtConstContour.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('ct_const.exe')
                self.plotCtConst(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(2), w_JSON_FILE, r_JSON_FILE, fname_exe)
            if window.radioButtonCt.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('ct.exe')
                self.plotCtKarman(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(1), w_JSON_FILE, r_JSON_FILE, fname_exe)
            if window.radioButtonCtContour.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('ct.exe')
                self.plotCtKarman(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(2), w_JSON_FILE, r_JSON_FILE, fname_exe)
            if window.radioButtonKarman.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('karman.exe')
                self.plotCtKarman(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(1), w_JSON_FILE, r_JSON_FILE, fname_exe)
            if window.radioButtonKarmanContour.isChecked():
                n, a, b, h, q0, E, nu, numIters, fname_exe = self.getValues('karman.exe')
                self.plotCtKarman(n, a, b, h, q0, E, nu, numIters, self.getXY(n, a, b), \
                    PlotType(2), w_JSON_FILE, r_JSON_FILE, fname_exe)

    def isChangedInputData(self, fname_exe):
        for i in range(len(self.inputData)):
            if self.getValues(fname_exe)[i] != self.inputData[i]:
                return True
        return False

    def startProgram(self, fname1, fname2, n, a, b, h, q0, E, nu, numIters, cmd, fname_exe):
        self.inputData = (n, a, b, h, q0, E, nu, numIters, fname_exe)
        self.deleteJSONFiles(fname1, fname2)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
        p.wait()
        return p

    def plotContours(self, ax, x, y, z, b):
        ax.contour(x, y, z, zdir='z', offset=0, cmap=cm.coolwarm, antialiased=True)
        ax.contour(x, y, z, zdir='x', offset=0, cmap=cm.coolwarm, antialiased=True)
        ax.contour(x, y, z, zdir='y', offset=b, cmap=cm.coolwarm, antialiased=True)

    def setAxes(self, ax, a, b, z):
        ax.set_xlabel('x')
        ax.set_xlim(0, a)
        ax.set_ylabel('y')
        ax.set_ylim(0, b)
        ax.set_zlabel(z)

    def plotSurfaceOrContour(self, fig, bool_const, plotType, a, b, x, y, z, z_str, num=0, text=''):
        if bool_const:
            ax = fig.gca(projection='3d')
        else:
            ax = fig.add_subplot(1, 2, num, projection='3d')
        self.setAxes(ax, a, b, z_str)
        if num == 1:
            ax.text2D(-0.2, -0.05, text, transform=ax.transAxes)
        if plotType.name == 'surface':
            ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=True)
        elif plotType.name == 'contour':
            self.plotContours(ax, x, y, z, b)

    def plotCtConst(self, n, a, b, h, q0, E, nu, numIters, tuple_xy, plotType, fname1, fname2, fname_exe):
        if os.path.isfile(fname2) or self.isChangedInputData(fname_exe) or not os.path.isfile(fname1):
            cmd = '{0} {1} {2} {3} {4} {5} {6} {7}'.format(fname_exe, n, a, b, h, q0, E, nu)
            self.startProgram(fname1, fname2, n, a, b, h, q0, E, nu, numIters, cmd, fname_exe)
        w = getListFromJSON(fname1, 'w')
        if w:
            z = np.array([[w[i][j] for i in range(n + 1)] for j in range(n + 1)])
            try:
                x = tuple_xy[0]
                y = tuple_xy[1]
                fig = plt.figure('Прогиб пластины')
                self.plotSurfaceOrContour(fig, True, plotType, a, b, x, y, w, 'w')
                plt.show()
            except (ValueError, TypeError):
                QMessageBox.warning(self, 'Ошибка','Поменяйте количество отрезков сетки')

    def plotCtKarman(self, n, a, b, h, q0, E, nu, numIters, tuple_xy, plotType, fname1, fname2, fname_exe):
        if not os.path.isfile(fname2) or self.isChangedInputData(fname_exe) or \
            (not os.path.isfile(fname1) and not os.path.isfile(fname2)):
            cmd = '{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(fname_exe, n, a, b, h, q0, E, nu, numIters)
            self.startProgram(fname1, fname2, n, a, b, h, q0, E, nu, numIters, cmd, fname_exe)
        _w = getListFromJSON(fname1, 'w')
        _r = getListFromJSON(fname2, 'r')
        if _w != None and _r != None: 
            w = np.array([[_w[i][j] for i in range(n + 1)] for j in range(n + 1)])
            r = np.array([[_r[i][j] for i in range(n + 1)] for j in range(n + 1)])
            try:
                x = tuple_xy[0]
                y = tuple_xy[1]
                fig = plt.figure(figsize=plt.figaspect(0.5), num='Прогиб и реакция пластины')
                if numIters == 0:
                    text = 'n = {0}\na = {1} см\nb = {2} см\nh = {3} см\nq0 = {4} кг/см^2\nE = {5} кг/см^2\n\u03bd = {6}' \
                        .format(n, a, b, h, q0, E, nu)
                else:
                    text = 'n = {0}\na = {1} см\nb = {2} см\nh = {3} см\nq0 = {4} кг/см^2\nE = {5} кг/см^2\n\u03bd = {6}\n\nколичество итераций = {7}' \
                        .format(n, a, b, h, q0, E, nu, numIters)
                self.plotSurfaceOrContour(fig, False, plotType, a, b, x, y, w, 'w', 1, text)
                self.plotSurfaceOrContour(fig, False, plotType, a, b, x, y, r, 'r', 2, text)
                plt.show()
            except (ValueError, TypeError):
                QMessageBox.warning(self, 'Ошибка','Поменяйте количество отрезков сетки')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())