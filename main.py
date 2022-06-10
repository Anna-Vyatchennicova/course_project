import sys
import os
import pandas as pd
import math
import cmath
import mpmath
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QMessageBox
from scipy import interpolate
from scipy import optimize
import numpy as np
from design import Ui_MainWindow
from functools import partial

m = n = R_x = F = x = I = delta_o = k_x = H_x = P_x = T_x = None


def func(x, R):
    a = mpmath.besselj(2, x * R)
    return float(a)


def f_for_v(sigm, ind):
    global R_x, F
    return 1 / math.pi * (I / (m * R_x[ind])) ** 2 * (F[1][ind] -
                                                    F[2][ind] * delta_o[ind] / T_x[ind]) ** 2 / (F[2][ind] * sigm ** 2)


def cubic_interp(A, B, x):
    tck = interpolate.splrep(A, B)
    return interpolate.splev(x, tck)


def find_in_interv(a, b, R):
    b += 10e-10
    try:
        while True:
            if a > b:
                raise ValueError
            sol = optimize.root_scalar(func, args=(R,), x0=a, x1=b, method='secant')
            if a <= sol.root <= b:
                if sol.converged:
                    return sol.root
                return None
            if sol.root > b:
                b -= 1
            else:
                a += 1
    except ValueError:
        return None


def f_tr(v, ind):
    global m, n, I, F, delta_o, k_x, H_x, P_x
    beta = R_x[ind] * cmath.sqrt(k_x[ind] / v)
    b = mpmath.besselj(2, beta) / mpmath.besselj(1, beta)
    i = nu = j = 0
    while i <= 20:
        nu = find_in_interv(nu + 300, nu + 700, R_x[ind])
        znak = cmath.sqrt(nu ** 2 - (beta / R_x[ind]) ** 2)
        y = cmath.tanh(znak * H_x[ind] / 2) / (nu ** 2 * znak ** 3)
        j += y
        i += 1
    sm = j
    rassch = -2 * v * m * b * beta + 4 * math.pi * n * P_x[ind] * R_x[ind] ** 2 * k_x[ind] ** 2 / v * sm
    f = rassch.real + F[1][ind] / (2 * math.pi) * rassch.imag - 2 * I * (F[1][ind] / F[2][ind] - delta_o[ind] / T_x[ind])
    return f


def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


class MyWidget(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.filename1 = self.filename2 = ''
        self.pushButton.clicked.connect(partial(self.choose_file, 1))
        self.pushButton_2.clicked.connect(partial(self.choose_file, 2))
        self.pushButton_3.clicked.connect(self.start)

    def start(self):
        global m, n, R_x, F, x, I, delta_o, k_x, H_x, P_x, T_x
        try:
            m = float(self.mass.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле массы образца не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            R = float(self.radius.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле радиуса образца не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            p = float(self.plotnost.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле плотности расплава не заполнено или заполнено некорректно!',
                                 QMessageBox.Ok)
            return
        try:
            I = float(self.moment_inerts.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле момента инерции полвесной системы'
                                                 ' не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            vmin = float(self.minim.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле минимального значения не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            vmax = float(self.max.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле максимального значения не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            t0 = float(self.period.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле периода пустой подвесной'
                                                 ' системы не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        try:
            delta_b = float(self.decrement.text().replace(',', '.'))
        except Exception:
            QMessageBox.critical(self, 'Ошибка', 'Поле декремента пустой подвесной '
                                                 'системы не заполнено или заполнено некорректно!', QMessageBox.Ok)
            return
        if self.torts_1.isChecked():
            n = 1
        else:
            n = 2
        if len(self.filename1) == 0:
            QMessageBox.critical(self, 'Ошибка', 'Не выбран файл с результатами измерений!', QMessageBox.Ok)
            return
        if len(self.filename2) == 0:
            QMessageBox.critical(self, 'Ошибка', 'Не выбран файл с данными для расчёта вязкости!', QMessageBox.Ok)
            return
        filename = self.filename1
        filename_s = self.filename2
        self.pushButton_3.setEnabled(False)
        try:
            if filename.endswith('xlsx'):
                file = pd.ExcelFile(filename, engine='openpyxl')
                sheet = file.sheet_names[0] sheet -
                stolb = pd.read_excel((filename), sheet_name=sheet, usecols='A', engine='openpyxl', header=None)
                amount_of_measurings = len(stolb)
                if isinstance(stolb[0][0], str):
                    flag = True
                    amount_of_measurings -= 1
                else:
                    flag = False
                data = pd.read_excel((filename), sheet_name=sheet, usecols='C:E', engine='openpyxl',
                                     header=(None if not flag else 0))
            else:
                file = pd.ExcelFile(filename, engine='xlrd')
                sheet = file.sheet_names[0]
                stolb = pd.read_excel((filename), sheet_name=sheet, usecols='A', engine='xlrd', header=None)
                amount_of_measurings = len(stolb)
                if isinstance(stolb[0][0], str):
                    flag = True
                    amount_of_measurings -= 1
                else:
                    flag = False
                data = pd.read_excel((filename), sheet_name=sheet, usecols='C:E', engine='xlrd',
                                     header=(None if not flag else 0))
            F = [np.array(data[key]) for key in data.keys()]
            usr = 10
            #alpha = 0.0000106
            x = range(0, amount_of_measurings)
            delta_o = [delta_b for e in x]
            T_x = [t0 for e in x]
            P_x = [p for e in x]
            R_x = [R for e in x]
            H_x = [m / (math.pi * R_x[i] ** 2 * P_x[i]) for i in x]
            vmin_x = [vmin for e in x]
            vmax_x = [vmax for e in x]
            k_x = [(F[1][i] + 2 * math.pi * 1j) / F[2][i] for i in x]
            f_min = [f_tr(vmin, i) for i in x]
            f_max = [f_tr(vmax, i) for i in x]
            v_x = []
            for index in x:
                if abs(f_max[index] + f_min[index]) < abs(f_max[index]) + abs(f_min[index]):
                    while ((vmax_x[index] - vmin_x[index]) / vmax_x[index]) > 0.0001:
                        vx = (vmax_x[index] + vmin_x[index]) / 2
                        err = f_tr(vx, index)
                        if err > 0:
                            vmax_x[index] = vx
                        elif err < 0:
                            vmin_x[index] = vx
                    v_x.append(vx)
                else:
                    v_x.append(0)
            if filename_s.endswith('xlsx'):
                file = pd.ExcelFile(filename_s, engine='openpyxl')
                sheet = file.sheet_names[0]
                stolb = pd.read_excel((filename_s), sheet_name=sheet, usecols='A', engine='openpyxl', header=None)
                amount_of_measurings = len(stolb)
                if isinstance(stolb[0][0], str):
                    flag = True
                    amount_of_measurings -= 1
                else:
                    flag = False
                data = pd.read_excel((filename_s), sheet_name=sheet, usecols='A:C', engine='openpyxl',
                                     header=(None if not flag else 0))
            else:
                file = pd.ExcelFile(filename_s, engine='xlrd')
                sheet = file.sheet_names[0]
                stolb = pd.read_excel((filename_s), sheet_name=sheet, usecols='A', engine='xlrd', header=None)
                amount_of_measurings = len(stolb)
                if isinstance(stolb[0][0], str):
                    flag = True
                    amount_of_measurings -= 1
                else:
                    flag = False
                data = pd.read_excel((filename_s), sheet_name=sheet, usecols='A:C', engine='xlrd',
                                     header=(None if not flag else 0))
            Y = np.array(data[data.keys()[0]])
            B = np.array(data[data.keys()[1]])
            C = np.array(data[data.keys()[2]])
            sigma_x = [1 for e in x]
            vh_x = []
            for ind in x:
                qq = 0
                vh = f_for_v(sigma_x[ind], ind)
                while (vh - qq) / vh > 0.0001:
                    y_x = R_x[ind] ** 2 * 2 * math.pi / (F[2][ind] * vh)
                    Xh_x = F[1][ind] / (2 * math.pi)
                    A_x = 3 / (math.sqrt(2 * y_x))
                    Bh_x = cubic_interp(Y, B, y_x)
                    Ch_x = cubic_interp(Y, C, y_x)
                    qq = vh
                    SIG = 1 - 1.5 * Xh_x - 3 / 8 * Xh_x ** 2 - A_x + 2 * n * R_x[ind] / H_x[ind] * (Bh_x - Ch_x * Xh_x)
                    vh = f_for_v(SIG, ind)
                y_x = R_x[ind] ** 2 * 2 * math.pi / (F[2][ind] * vh)
                if y_x < 100:
                    vh = 0
                vh_x.append(vh)
            Rez = []
            for e in F:
                Rez.append(e)
            Rez.append(np.array(v_x))
            Rez.append(np.array(vh_x))
            pp = int((len(x) + 1) / usr - 1)
            Rez_sr = []
            for pos in range(5):
                stolb = []
                for ind in range(0, pp):
                    j = 0
                    i = ind * usr
                    while ind * usr <= i < (ind + 1) * usr:
                        j += Rez[pos][i]
                        i += 1
                    stolb.append(j / usr)
                Rez_sr.append(stolb[::])
            df = pd.DataFrame({i: Rez[i] for i in range(5)})
            df.to_excel('./data2.xlsx', index=False)
            df = pd.DataFrame({i: Rez_sr[i] for i in range(5)})
            df.to_excel('./data2sr.xlsx', index=False)
            QMessageBox.information(self, 'Успех', 'Расчёт прошёл успешно. Таблицы с '
                                                   'результатами доступны в папке с программой', QMessageBox.Ok)
        except Exception as e:
            print(e)
            QMessageBox.critical(self, 'Ошибка', 'Что-то пошло не так... Проверьте корректоность данных',
                                 QMessageBox.Ok)
        self.pushButton_3.setEnabled(True)

    def choose_file(self, a0):
        passwords_file = QFileDialog.getOpenFileName(self, 'Выбор файла', '',
                                                     "Excel files (*.xls; *.xlsx)")[0]
        if len(passwords_file) > 0:
            if a0 == 1:
                self.filename1 = passwords_file
                test_name = '../' + '/'.join(passwords_file.split('/')[-2:])
                if len(test_name) <= 24:
                    self.pushButton.setText(test_name)
                else:
                    self.pushButton.setText('../' + test_name.split('/')[-1][:17])
            else:
                self.filename2 = passwords_file
                test_name = '../' + '/'.join(passwords_file.split('/')[-2:])
                if len(test_name) <= 24:
                    self.pushButton_2.setText(test_name)
                else:
                    print(len('../' + test_name.split('/')[-1][:22]))
                    self.pushButton_2.setText('../' + test_name.split('/')[-1][:17])


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyWidget()
    ex.setFixedSize(930, 635)
    ex.setWindowIcon(QIcon(resource_path('python_pythonlover.ico')))
    ex.show()
    sys.exit(app.exec_())


#stolbi = pd.read_excel((self.skl), sheet_name='Интеграл_спектры ТЗЧ СКЛ и ГКЛ', usecols='C', header=6)