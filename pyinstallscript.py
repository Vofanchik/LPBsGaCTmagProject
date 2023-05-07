import PyInstaller.__main__

PyInstaller.__main__.run([
    'main.py',
    '--onefile',
    '--windowed',
    '-iexternal/resources/Img/breast-cancer_cell-transformed.ico'
])