import os

import PyInstaller.__main__
os.environ["R_HOME"] = 'C:\Program Files\R\R-4.3.1'
PyInstaller.__main__.run([
    'main.py',
    '--onefile',
    '--windowed',
    '-iexternal/resources/Img/breast-cancer_cell-transformed.ico'
])