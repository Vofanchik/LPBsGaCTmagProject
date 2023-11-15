import base64
from datetime import datetime
from PyQt6 import QtWidgets
from PyQt6.QtCore import QByteArray
from PyQt6.QtGui import QAction, QMovie, QIcon, QPixmap, QImage
from PyQt6.QtSvg import QSvgRenderer
from PyQt6.QtSvgWidgets import QSvgWidget, QGraphicsSvgItem
from PyQt6.QtWidgets import QDialog, QTableWidgetItem, QSizePolicy, QLabel, QTextEdit, QMenu, QDateEdit, \
    QFileDialog, QComboBox, QMessageBox, QWidget, QGraphicsView, QGraphicsScene
from PyQt6.QtCore import Qt
from PyQt6.uic.properties import QtCore

from Classes.DataBase import DataBase
from Classes.Testin import timeit
from Classes.dateclass import from_dot_to_rec
from Classes.JchemPaint import runJCP
from Classes.RDkit import getMolSvg, similiaryty_list_return, similiaryty_map_rerurn, similiaryty_map_rerurn_png
from Classes.Rclass import calculate_EC50_SE_plots

from UI_files.DialogAddNewCompound import Ui_Dialog
from UI_files.MainWindow import Ui_MainWindow
from UI_files.MTTable import Ui_MTTable




class SvgShow(QDialog):
    def __init__(self, svg):
        super(SvgShow, self).__init__()

        self.setMinimumWidth(400)
        self.setMaximumHeight(400)
        self.layout = QtWidgets.QGridLayout()
        self.setWindowTitle("Карта сходства соединений")
        self.setLayout(self.layout)

        policy = QSizePolicy()
        policy.setWidthForHeight(True)

        self.svg_widget = QSvgWidget()
        self.svg_widget.setSizePolicy(policy)

        self.svg_widget.setMinimumSize(800, 800)
        self.svg_widget.load(QByteArray(bytes(svg, encoding='utf-8')))




        self.layout.addWidget(self.svg_widget)


class ShowPng(QDialog):
    def __init__(self, png):
        super().__init__()

        # self.setGeometry(200, 200, 200, 200)

        label = QLabel(self)
        # label.setMaximumHeight(100)
        qp = QPixmap()
        qp.loadFromData(png)
        label.setPixmap(qp)

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(label)
        self.setWindowTitle("Карта сходства соединений")
        self.setLayout(self.layout)



class MTTtableWidget(QDialog):
    def __init__(self, mode=False, **kwargs):
        super(MTTtableWidget, self).__init__()
        self.kwargs = kwargs
        self.ui = Ui_MTTable()
        self.ui.setupUi(self)
        self.db = DataBase()

        if not mode:
            self.fill_MTT_table()
        elif mode=='Similarity':
            self.init_similarity_mode()


    def fill_MTT_table(self):
        for co, it in enumerate(self.db.show_mtt(self.kwargs['chosen_compound_id'])):
            self.ui.tableWidget.setRowCount(co + 1)
            self.ui.tableWidget.setItem(co, 0, QTableWidgetItem(it[0]))
            self.ui.tableWidget.setItem(co, 1, QTableWidgetItem(from_dot_to_rec(1,it[1])))
            self.ui.tableWidget.setItem(co, 3, QTableWidgetItem(it[14]))
            self.ui.tableWidget.setItem(co, 2, QTableWidgetItem(it[7]))

            self.ui.tableWidget.setItem(co, 4, QTableWidgetItem(str(it[4])))
            self.ui.tableWidget.setItem(co, 5, QTableWidgetItem(str(it[5])))

    def init_similarity_mode(self):
        self.setWindowTitle(f"Подобия {self.db.show_name_by_id(self.kwargs['chosen_compound_id'])[0]}")
        [self.ui.tableWidget.setColumnHidden(i, True) for i in range(4,6)]
        self.ui.tableWidget.setColumnHidden(0, True)
        [self.ui.tableWidget.horizontalHeaderItem(n+1).setText(i) for n,i in enumerate(
            ['Наименование', '% сходства', 'EC50'])]
        self.ui.tableWidget.setColumnWidth(1, 400)
        self.fill_similiarity_table()
        self.ui.tableWidget.itemDoubleClicked.connect(self.show_similiarity_map)

    def fill_similiarity_table(self):
        for co, it in enumerate(similiaryty_list_return(self.db.show_smile_by_id(self.kwargs['chosen_compound_id'])[0],
                                                        self.db.show_all_id_name_smile())[0:10]):
            self.ui.tableWidget.setRowCount(co + 1)
            self.ui.tableWidget.setItem(co, 0, QTableWidgetItem(str(it[0])))
            self.ui.tableWidget.setItem(co, 1, QTableWidgetItem(it[1]))
            self.ui.tableWidget.setItem(co, 2, QTableWidgetItem(str(round(it[2]*100,2))))
            ec=str(self.db.show_best_mtt_result(it[0])[0])
            self.ui.tableWidget.setItem(co, 3, QTableWidgetItem(ec if ec !="None" else ""))


    def show_similiarity_map(self):
        targ = self.db.show_smile_by_id(self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text())
        ref = self.db.show_smile_by_id(self.kwargs['chosen_compound_id'])


        # dlg = SvgShow(similiaryty_map_rerurn(targ[0], ref[0]))
        dlg = ShowPng(similiaryty_map_rerurn_png(targ[0], ref[0]))

        dlg.exec()


class DialogAddNewCompound(QDialog):
    def __init__(self, mode=False, **kwargs):
        self.kwargs = kwargs
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
        elif  mode == "add line":
            self.add_line_mode()
        elif  mode == "add MTT":
            self.add_mtt_mode()

    def getNameFromArea(self):
        return self.ui.lineEdit.text()

    def getSheetFromCombo(self):
        return self.combo_sheet.currentText()

    def getDateFromCalendar(self):
        return self.calendar_widget.text()

    def getLineNameFromCombo(self):
        return self.ui.comboBox.currentText()

    def getDescriptFromTextEdit(self):
        return self.line_description.toPlainText()

    def getSolFromCombo(self):
        ls = list(i[1] for i in self.db.show_solutors())
        i = ls.index(self.ui.comboBox.currentText())
        return self.db.show_solutors()[i][0]

    def getLineFromCombo(self):
        ls = list(i[1] for i in self.db.show_lines())
        i = ls.index(self.ui.comboBox.currentText())
        return self.db.show_lines()[i][0]

    def add_sol_mode(self):
        self.setWindowTitle("Добавить растворитель")
        self.ui.comboBox.hide()
        self.ui.label_2.hide()

    def add_line_mode(self):
        self.setWindowTitle("Добавить линию")
        self.ui.comboBox.hide()
        self.ui.label_2.hide()
        self.ui.gridLayout.addWidget(QLabel(text="Описание"), 0, 3, alignment=Qt.AlignmentFlag.AlignCenter)
        self.line_description = QTextEdit()
        self.ui.gridLayout.addWidget(self.line_description, 1, 3, alignment=Qt.AlignmentFlag.AlignCenter)

    def insert_filepath(self):
        self.file_path = QFileDialog.getOpenFileName(None, "Select file", "", "xlsx files (*.xlsx)")
        self.ui.lineEdit.setText(self.file_path[0])
        self.ui.lineEdit.setDisabled(True)

    def add_mtt_mode(self):
        self.setWindowTitle("Добавить МТТ")
        self.ui.label.setText("Путь к файлу")
        self.ui.label_2.setText("Клеточная линия")
        self.ui.lineEdit.mouseDoubleClickEvent = lambda event: self.insert_filepath()
        self.ui.comboBox.addItems(i[1] for i in self.db.show_lines())
        self.calendar_widget = QDateEdit()
        self.combo_sheet = QComboBox()
        self.combo_sheet.addItems([str(i) for i in range(1,6)])
        self.calendar_widget.setCalendarPopup(True)
        self.calendar_widget.setDate(datetime.now())
        self.ui.gridLayout.addWidget(QLabel('Номер листа книги xls'), 4,2)
        self.ui.gridLayout.addWidget(self.combo_sheet, 5, 2)
        self.ui.gridLayout.addWidget(QLabel('Дата проведения'), 6,2)
        self.ui.gridLayout.addWidget(self.calendar_widget, 7,2)
        self.ui.gridLayout.addWidget(self.ui.buttonBox, 8, 2)
        self.ui.gridLayout.addWidget(QLabel(self.kwargs['chosen_compound']), 9, 2)
        # self.ui.gridLayout.removeWidget(self.ui.buttonBox)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.id = None
        self.db = DataBase()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.load_dialog = QDialog()
        self.load_dialog_layout = QtWidgets.QVBoxLayout()
        self.load_dialog.setWindowTitle('Подождите')
        self.load_dialog.setGeometry(500,200,1,1)
        self.load_dialog.setLayout(self.load_dialog_layout)
        self.load_dialog_label = QtWidgets.QLabel()
        self.load_dialog_movie = QMovie("./external/resources/gif/load.gif")
        self.load_dialog_label.setMovie(self.load_dialog_movie)
        self.load_dialog_layout.addWidget(self.load_dialog_label)

        try:
            self.setStyleSheet(open("./Ui_files/Style.qss", "r").read())
        except:
            pass
        self.ui.tableWidget.setColumnHidden(0, 1)
        # self.ui.tableWidget.setColumnHidden(5, 1)
        self.ui.pushButton.clicked.connect(self.add_compound)
        self.ui.delCompoundButton.clicked.connect(self.del_compound)
        self.fill_table_compounds()

        policy = QSizePolicy()
        policy.setWidthForHeight(True)

        self.svg_widget = QSvgWidget()
        self.svg_widget.setSizePolicy(policy)

        self.svg_widget.setMinimumSize(300, 300)

        def add_menu():
            bar = self.menuBar()
            actions = bar.addMenu("Действия")
            actions.addAction("Добавить растворитель")
            actions.addAction("Добавить линию")
            # open_groups = QAction("Open groups", self)
            actions.triggered[QAction].connect(self.menu_bar_triggered)

        add_menu()
        self.ui.tableWidget.itemDoubleClicked.connect(self.enter_chosen_compound)
        self.contextMenu = QMenu(self)
        self.addMTTAction = QAction("Добавить МТТ-тест", self)
        self.ShowMTTByIdAction = QAction("Показать тесты", self)
        self.ShowSimiliarity = QAction("Найти подобия", self)
        self.contextMenu.addAction(self.addMTTAction)
        self.contextMenu.addAction(self.ShowMTTByIdAction)
        self.contextMenu.addAction(self.ShowSimiliarity)
        self.addMTTAction.setDisabled(True)
        self.ShowMTTByIdAction.setDisabled(True)
        self.ShowSimiliarity.setDisabled(True)
        self.ui.tableWidget.cellClicked.connect(lambda: [self.ui.delCompoundButton.setEnabled(True),
                                                         self.addMTTAction.setEnabled(True),
                                                         self.ShowMTTByIdAction.setEnabled(True),
                                                         self.ShowSimiliarity.setEnabled(True)])

        self.addMTTAction.triggered.connect(self.add_mtt_dialog)
        self.ShowMTTByIdAction.triggered.connect(self.show_mtt_dialog)
        self.ShowSimiliarity.triggered.connect(self.show_similiarity_dialog)
        self.ui.tableWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.ui.tableWidget.customContextMenuRequested.connect(self.showContextMenu)

    def showContextMenu(self, position):
        self.contextMenu.exec(self.ui.tableWidget.mapToGlobal(position))

    def add_mtt_dialog(self):
        chosen_compound = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 1).text()
        chosen_compound_id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        dlg = DialogAddNewCompound(mode="add MTT",
                                   chosen_compound=chosen_compound,
                                   chosen_compound_id = chosen_compound_id)
        file_path = ''
        date_test = ''
        sheet_num = ''
        line_num = ''
        line_name = ''
        rej = False

        if dlg.exec():
            rej = dlg.rejected
            file_path = dlg.getNameFromArea()
            date_test = dlg.getDateFromCalendar()
            sheet_num = dlg.getSheetFromCombo()
            line_num = dlg.getLineFromCombo()
            line_name = dlg.getLineNameFromCombo()

        if rej is False:
            return



        self.load_dialog.show()
        self.startAnimation()
        try:
            ic_ec = calculate_EC50_SE_plots(file_path,chosen_compound, line_name, int(sheet_num))
            self.db.add_MTT(from_dot_to_rec(0, date_test), chosen_compound_id, line_num, round(ic_ec[0], 2),
                            round(ic_ec[1], 2))
            self.load_dialog.hide()

        except:
            self.load_dialog.hide()
            QMessageBox().critical(self, 'Ошибка', "Что-то пошло не так")

    def show_mtt_dialog(self):
        chosen_compound_id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        dlg = MTTtableWidget(chosen_compound_id=chosen_compound_id)
        dlg.ui.tableWidget.setColumnHidden(0, True)

        if dlg.exec():
            pass

    def show_similiarity_dialog(self):
        chosen_compound_id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        dlg = MTTtableWidget(mode='Similarity', chosen_compound_id=chosen_compound_id)

        dlg.exec()




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

        elif press.text() == "Добавить линию":
            dlg = DialogAddNewCompound(mode="add line")
            val = ''
            descript = ''
            rej = False

            if dlg.exec():
                rej = dlg.rejected
                val = dlg.getNameFromArea()
                descript = dlg.getDescriptFromTextEdit()
            print(descript)

            if rej is False:
                return

            self.db.add_cellLine(val, descript)

    def startAnimation(self):
        self.load_dialog_movie.start()

    def stopAnimation(self):
        self.load_dialog_movie.stop()

    def enter_chosen_compound(self):
        self.id = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 0).text()
        self.smiles = self.ui.tableWidget.item(self.ui.tableWidget.currentRow(), 2).text()

        try:
            self.ui.gridLayout.removeWidget(self.ui.gridLayout.itemAtPosition(1,2).widget())
        except:
            pass


        self.ui.gridLayout.addWidget(self.svg_widget, 0, 2, alignment=Qt.AlignmentFlag.AlignCenter)
        # self.ui.gridLayout.addWidget(QtWidgets.QLabel
        #                              (self.ui.tableWidget.item
        #                               (self.ui.tableWidget.currentRow(), 5).text()), 1, 2)
        self.update_svg_widget(self.smiles)

    @timeit
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
        self.db.add_compound(name=name, smiles=smile, solutor_id=solutor)
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
            self.ui.tableWidget.setItem(co, 4, QTableWidgetItem(f"{it[8]}"))
            self.ui.tableWidget.setItem(co, 5, QTableWidgetItem(f"{it[5]}"))
            self.ui.tableWidget.setItem(co, 3, QTableWidgetItem(str(round(it[4], 2)) if type(it[4])==float else str(it[4])))


