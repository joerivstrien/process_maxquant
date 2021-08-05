import PyInstaller.__main__
import os

workdir = os.getcwd()
dist_path = os.path.join(workdir, 'dist')
print(dist_path)
PyInstaller.__main__.run([
    'gui_file_acceptor.py',
    '--distpath', dist_path,
    '--onefile',
    '--name', "process_maxquant_executable_with_gui_3"])
