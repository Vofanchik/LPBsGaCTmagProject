import sys

from PyQt6.QtGui import QIcon, QPixmap
from PyQt6.QtWidgets import QApplication
from Classes.MainWindow import MainWindow

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowIcon(QIcon(QPixmap('external/resources/Img/breast-cancer_cell-transformed.png')))
    mw.show()

    sys.exit(app.exec())