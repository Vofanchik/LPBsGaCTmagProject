from PyQt6 import QtWidgets
from PyQt6.QtCore import QByteArray, Qt
from PyQt6.QtGui import QAction
from PyQt6.QtSvgWidgets import QSvgWidget
from PyQt6.QtWidgets import QDialog, QTableWidgetItem, QSizePolicy
from PyQt6.QtCore import Qt

from Classes.DataBase import DataBase
from Classes.JchemPaint import runJCP
from Classes.RDkit import getMolSvg, iupac_from_smiles
from UI_files.DialogAddNewCompound import Ui_Dialog
from UI_files.MainWindow import Ui_MainWindow


class DialogAddNewCompound(QDialog):
    def __init__(self, mode=False):
        super(DialogAddNewCompound, self).__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.db = DataBase()

        if not mode:
            self.ui.comboBox.addItems(i[1] for i in self.db.show_solutors())
            self.setWindowTitle("Добавить соединение")

        try:
            self.setStyleSheet(open("./Ui_files/Style.qss", "r").read())
        except:
            pass

        if mode == "add solutor":
            self.add_sol_mode()

    def getNameFromArea(self):
        return self.ui.lineEdit.text()

    def getSolFromCombo(self):
        ls = list(i[1] for i in self.db.show_solutors())
        i = ls.index(self.ui.comboBox.currentText())
        return self.db.show_solutors()[i][0]

    def add_sol_mode(self):
        self.setWindowTitle("Добавить растворитель")
        self.ui.comboBox.hide()
        self.ui.label_2.hide()


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.id = None
        self.db = DataBase()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        try:
            self.setStyleSheet(open("./Ui_files/Style.qss", "r").read())
        except:
            pass
        self.ui.tableWidget.setColumnHidden(0, 1)
        self.ui.tableWidget.setColumnHidden(5, 1)
        self.ui.pushButton.clicked.connect(self.add_compound)
        self.ui.delCompoundButton.clicked.connect(self.del_compound)
        self.fill_table_compounds()

        policy = QSizePolicy()
        policy.setWidthForHeight(True)

        self.svg_widget = QSvgWidget()
        self.svg_widget.setSizePolicy(policy)


        self.svg_widget.setMinimumSize(300, 300)

        self.ui.tableWidget.cellClicked.connect(lambda: self.ui.delCompoundButton.setEnabled(True))

        # self.ui.gridLayout.setColumnMinimumWidth(1, 300)
        # self.ui.gridLayout.setColumnStretch(1,9)
        # self.ui.gridLayout.setRowMinimumHeight(0, 300)
        # self.ui.gridLayout.setRowStretch(0,9)

        def add_menu():
            bar = self.menuBar()
            actions = bar.addMenu("Действия")
            actions.addAction("Добавить растворитель")

            # open_groups = QAction("Open groups", self)
            actions.triggered[QAction].connect(self.menu_bar_triggered)

        add_menu()

        self.ui.tableWidget.itemDoubleClicked.connect(self.enter_chosen_compound)

    def menu_bar_triggered(self, press):
        if press.text() == "Добавить растворитель":
            dlg = DialogAddNewCompound(mode="add solutor")
            val = ''
            rej = False

            if dlg.exec():
                rej = dlg.rejected
                val = dlg.getNameFromArea()

            if rej is False:
                return

            self.db.add_solutor(val)

    def enter_chosen_compound(self):
        self.id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        self.smiles = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 2).text()

        try:
            self.ui.gridLayout.removeWidget(self.ui.gridLayout.itemAtPosition(1,2).widget())
        except:
            pass


        self.ui.gridLayout.addWidget(self.svg_widget, 0, 2, alignment=Qt.AlignmentFlag.AlignCenter)
        self.ui.gridLayout.addWidget(QtWidgets.QLabel
                                     (self.ui.tableWidget.item
                                      (self.ui.tableWidget.currentRow(), 5).text()), 1, 2)
        self.update_svg_widget(self.smiles)


    def update_svg_widget(self, smile):
        self.svg_widget.load(QByteArray(bytes(getMolSvg(smile), encoding='utf-8')))

    def add_compound(self):
        runJCP()
        dlg = DialogAddNewCompound()
        name = ''
        solutor = 0
        rej = False


        if dlg.exec():
            rej = dlg.rejected
            name = dlg.getNameFromArea()
            solutor = dlg.getSolFromCombo()

        if rej is False:
            return

        f = open('./smile.txt', 'r')
        smile = f.read()
        f.close()
        self.db.add_compound(name, smile, solutor)
        self.fill_table_compounds()
        self.ui.statusbar.showMessage("Соединение добавлено", 5000)

    def del_compound(self):
        self.id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        self.db.del_compound_by_id(self.id)
        self.fill_table_compounds()
        self.ui.statusbar.showMessage("Соединение удалено", 5000)

    def fill_table_compounds(self):
        for co, it in enumerate(self.db.show_compounds()):
            self.ui.tableWidget.setRowCount(co + 1)
            self.ui.tableWidget.setItem(co, 0, QTableWidgetItem(f"{it[0]}"))
            self.ui.tableWidget.setItem(co, 1, QTableWidgetItem(f"{it[1]}"))
            self.ui.tableWidget.setItem(co, 2, QTableWidgetItem(f"{it[2]}"))
            self.ui.tableWidget.setItem(co, 4, QTableWidgetItem(f"{it[7]}"))
            self.ui.tableWidget.setItem(co, 5, QTableWidgetItem(f"{it[5]}"))
            self.ui.tableWidget.setItem(co, 3, QTableWidgetItem(f"{it[4]}"))

