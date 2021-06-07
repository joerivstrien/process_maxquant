import PyInstaller.__main__
import os

workdir = os.getcwd()
dist_path = os.path.join(workdir, 'dist')
print(dist_path)
PyInstaller.__main__.run([
    'process_maxquant.py',
    '--distpath', dist_path,
    '--onefile',
    '--name', "test_6"])
