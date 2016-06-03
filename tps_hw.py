#
#   Istanbul Technical University
#
#   Theory of Plates and Shells - UUM608E
#   Homework 4
#
#   Süleyman Muti
#   511142128
#
#   22.05.2016
#

import sys

import numpy as np
import matplotlib.pyplot as plt

from PyQt4 import QtGui
from PyQt4.uic import loadUiType

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


Ui_MainWindow, QMainWindow = loadUiType('tps_hw_GUI.ui')


class Main(QMainWindow, Ui_MainWindow):

    # Define plate arrays as global variables to be able use them elsewhere in the program.
    global w_plate
    global mx_plate
    global my_plate
    global qx_plate
    global qy_plate

    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)

        self.btnPlot.setStyleSheet('color: red')
        self.textEdit.setTextColor(QtGui.QColor("magenta"))

        # Check deflections radio button by default
        self.deflections_radio.setChecked(False)

        # Creating a QButtonGroup and adding the plot options ratio buttons.
        # This will be used to make sure new results are plotted every time the solve button is pressed.
        # Otherwise, if the deflections radio button is already selected, new run won't plot results.
        self.group = QtGui.QButtonGroup()
        self.group.addButton(self.deflections_radio)
        self.group.addButton(self.mx_radio)
        self.group.addButton(self.my_radio)
        self.group.addButton(self.qx_radio)
        self.group.addButton(self.qy_radio)
        self.plot_options_group.setEnabled(False)

        self.infoLabel.setStyleSheet('color: blue')
        self.infoLabel.setText('  Theory of Plates and Shells\n' +
                               '  Süleyman Muti\n' +
                               '  Spring 2016')

        self.btnPlot.clicked.connect(lambda: self.fdm())
        self.deflections_radio.toggled.connect(lambda: self.deflections_radio_checked())
        self.mx_radio.toggled.connect(lambda: self.mx_radio_checked())
        self.my_radio.toggled.connect(lambda: self.my_radio_checked())
        self.qx_radio.toggled.connect(lambda: self.qx_radio_checked())
        self.qy_radio.toggled.connect(lambda: self.qy_radio_checked())

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
                                         self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)

    def rmmpl(self, ):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def get_grid_size(self):
        if self.rBtn15x10.isChecked():
            return 15, 1e-6
        if self.rBtn30x20.isChecked():
            return 30, 1e-6
        elif self.rBtn180x120.isChecked():
            return 180, 1e-7
        elif self.rBtn540x360.isChecked():
            return 540, 1e-9

    def deflections_radio_checked(self):
        if self.deflections_radio.isChecked():
            self.rmmpl()
            self.print_plot(w_plate, 'Plate Deflections', 'viridis')

    def mx_radio_checked(self):
        if self.mx_radio.isChecked():
            self.rmmpl()
            self.print_plot(mx_plate, 'Mx Bending Moment', 'plasma')

    def my_radio_checked(self):
        if self.my_radio.isChecked():
            self.rmmpl()
            self.print_plot(my_plate, 'My Bending Moment', 'plasma')

    def qx_radio_checked(self):
        if self.qx_radio.isChecked():
            self.rmmpl()
            self.print_plot(qx_plate, 'Qx Transverse Shear Force', 'inferno')

    def qy_radio_checked(self):
        if self.qy_radio.isChecked():
            self.rmmpl()
            self.print_plot(qy_plate, 'Qy Transverse Shear Force', 'inferno')

    def print_plot(self, array_to_be_plotted, plot_title, colormap):
        a = array_to_be_plotted
        self.rmmpl()
        fig = Figure()
        self.addmpl(fig)
        ax = fig.add_subplot(111)
        cax = ax.imshow(a, cmap=colormap)
        fig.colorbar(cax, shrink=0.6, aspect=5)
        fig.tight_layout()
        ax.set_title(plot_title, fontsize=20, y=1.01)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.format_coord = lambda x, y: ''
        ax.set_xlabel('300 mm', fontsize=20, rotation=0)
        ax.set_ylabel('200 mm', fontsize=20, rotation=90)
        ax.autoscale(False)
        plt.show()

    def print_info(self, grid_size):
        grid_string = 'Grid size: {:d} x {:d}'.format(grid_size, grid_size*2//3)
        deflection_string = 'Maximum deflection of the plate: {0:.4f} mm'.format(np.max(w_plate))
        moment_string = 'Max. Mx: {:.4f} Nmm/mm \t Max. My: {:.4f} Nmm/mm'.format(np.max(mx_plate), np.max(my_plate))
        shear_string = 'Max. Qx:   {:.4f} N/mm \t Max. Qy:   {:.4f} N/mm'.format(np.max(qx_plate),
                                                                                        np.max(qy_plate))

        self.textEdit.setText(grid_string + '\n' + deflection_string + '\n' + moment_string + '\n' + shear_string)

    def fdm(self):

        # Let the function fdm() know that the plate arrays are global variables and can be used outside its scope.
        global w_plate
        global mx_plate
        global my_plate
        global qx_plate
        global qy_plate

        # Enabling initially disabled plot options area which was disabled to prevent
        # plot attempts without array creations.
        self.plot_options_group.setEnabled(True)

        # Deselecting the deflections radio button to make sure new results always get plotted.
        self.group.setExclusive(False)
        self.deflections_radio.setChecked(False)
        self.group.setExclusive(True)

        dim_x = 300  # Plate length in x-direction [mm]
        # dim_y = 200  # Plate length in y-direction [mm]

        # Use Python tuple to get the number of elements in the x-direction and corresponding convergence criterion.
        (m, conv) = self.get_grid_size()
        n = m * 2 // 3  # Number of elements in the y-direction.

        # Element size
        delta = dim_x / m

        # Initialize matrices for iterative solution.
        w_old = np.zeros((m + 3, n + 3))  # Initialize plate deflections to zero.
        w = np.copy(w_old)
        m1 = np.zeros((m + 3, n + 3))
        m2 = np.zeros((m + 3, n + 3))
        qx = np.zeros((m + 3, n + 3))
        qy = np.zeros((m + 3, n + 3))

        # Material properties, loading, and other properties
        e = 70000  # Modulus of elasticity [N/mm2]
        nu = 0.3  # Poisson's ratio
        h = 2  # Plate thickness [mm]
        po = 0.01  # Distributed load [N/mm2]

        # Plate stiffness
        d = (e * (h ** 3)) / (12 * (1.0 - (nu ** 2)))

        # Distributed load
        p = (po * (delta ** 4)) / d

        # Set logical condition to check if convergence criterion is met.
        cond = False  # Set to false to initiate iteration loop.

        # Loop to iterate plate deflections using Finite Difference Method

        iteration_number = 0  # Keep an eye on the iteration number.
        while not cond:
            cond = True  # set to true to check if criterion is met for all nodes.

            # Apply boundary conditions. Simply supported on all edges.
            w[:, 0] = -w_old[:, 2]
            w[:, n + 2] = -w_old[:, n]
            w[0, :] = -w_old[2, :]
            w[m + 2, :] = -w_old[m, :]

            # Calculate deflection of each node using neighbouring nodes
            # (i.e., using FDM)
            # Python range() statement is seemingly upper-bound exclusive. We need to add 1 to cover all range.
            for i in range(2, m + 1):

                # Check if the current index point has distributed load acting on it. If so apply the load, if not let it be zero.
                if i <= (m / 2 + 1):
                    k = p
                else:
                    k = 0

                for j in range(2, n + 1):

                    # Finite Difference Method Formula: v4 = p/d
                    w[i, j] = (1 / 20) * (k + 8 * (w[i + 1, j] + w[i - 1, j] + w[i, j + 1] + w[i, j - 1]) - 2 * (
                        w[i + 1, j + 1] + w[i - 1, j + 1] + w[i + 1, j - 1] + w[i - 1, j - 1]) - w[i + 2, j] - w[
                                              i - 2, j] - w[
                                              i, j + 2] - w[i, j - 2])
                    # Check if convergence criterion is met for each node. Set logical condition
                    # to false even if a single node violates the the condition so that the
                    # iteration can continue until all nodes meet the criterion.
                    if abs(w[i, j] - w_old[i, j]) > conv:
                        cond = False
            # Reset deflection matrices for next iteration.
            w_old = np.copy(w)

            # Keep an eye on the iteration number.
            iteration_number += 1

        # Calculate the bending moments and transverse shear forces based on the deflections.
        for i in range(2, m + 1):
            for j in range(2, n + 1):
                # Common terms in bending moment equations.
                m1[i, j] = -(w[i+1, j] - 2*w[i, j] + w[i-1, j])*d/delta**2
                m2[i, j] = -(w[i, j+1] - 2*w[i, j] + w[i, j-1])*d/delta**2

                # Transverse shear forces.
                qx[i, j] = -((w[i+2, j] - 2*w[i+1, j] + 2*w[i-1, j] - w[i-2, j]) +
                             (w[i+1, j+1] - 2*w[i+1, j] + w[i+1, j-1] - w[i-1, j+1] +
                              2*w[i-1, j] - w[i-1, j-1])) * d / (2*delta ** 3)
                qy[i, j] = -((w[i+1, j+1] - 2*w[i, j+1] + w[i-1, j+1] - w[i+1, j-1] + 2*w[i, j-1] - w[i-1, j-1]) +
                             (w[i, j+2] - 2*w[i, j+1] + 2*w[i, j-1] - w[i, j-2])) * d / (2*delta ** 3)

        # Assemble bending moment arrays.
        mx = m1 + nu*m2
        my = m2 + nu*m1

        # Exclude the ghost nodes that were necessary to apply boundary conditions,
        # and obtain deflections for plate nodes only. Plate deflections will be plotted correcting the orientation.
        w_plate = w[1:m + 2, 1:n + 2].transpose()
        mx_plate = mx[1:m + 2, 1:n + 2].transpose()
        my_plate = my[1:m + 2, 1:n + 2].transpose()
        qx_plate = qx[1:m + 2, 1:n + 2].transpose()
        qy_plate = qy[1:m + 2, 1:n + 2].transpose()

        # Set deflections radio button to checked to display the new run's results.
        self.deflections_radio.setChecked(True)

        # Print information summarizing the  solution
        # Maximum deflection, maximum bending moments, and maximum shear forces
        self.print_info(m)

        # Print results to console
        # print('\nW_max: ' + str(np.max(w_plate)))
        # print('Mx_max: ' + str(np.max(mx_plate)) + '\t' + 'My_max: ' + str(np.max(my_plate)))
        # print('Qx_max: ' + str(np.max(qx_plate)) + '\t' + 'Qy_max: ' + str(np.max(qy_plate)))

if __name__ == '__main__':
    fig1 = Figure()

    app = QtGui.QApplication(sys.argv)
    app.setStyle("plastique")
    main = Main()
    main.show()

    main.addmpl(fig1)

    sys.exit(app.exec_())
