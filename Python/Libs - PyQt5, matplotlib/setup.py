# -*- coding: UTF-8 -*-
'''
Модуль для сборки интерпретатора python и всех необходимых библиотек.
Для сборки использовать команду в терминале python setup.py build или bdist_msi 
'''

import sys, os
from cx_Freeze import setup, Executable

PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

setup(
    name = 'DeflectionRectPlate',
    version = '0.1',
    description = '''Программа построения поверхностей функций прогиба 
                и реакции прямоугольной пластины по различным теориям''',
    options = {'build_exe': {
        'includes': [
            'numpy.core._methods', 'numpy.lib.format', 'matplotlib.backends.backend_tkagg',
            'tkinter', 'tkinter.filedialog'
        ],
        'include_files':[
            os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
            os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll'),
            #r'C:\Windows\SysWOW64\api-ms-win-crt-stdio-l1-1-0.dll',
            'ct.exe',
            'ct_const.exe',
            'window.ui'
        ],
    }},
    executables = [
        Executable(script='DeflectionRectPlate.pyw', targetName='DeflectionRectPlate.exe', base=base)
    ]
)