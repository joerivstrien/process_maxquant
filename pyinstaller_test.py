import PyInstaller.__main__
import os

workdir = os.getcwd()
spec_path = os.path.join(workdir, 'process_maxquant.spec')
dist_path = os.path.join(workdir, 'dist')
print(spec_path, dist_path)
PyInstaller.__main__.run([
    'process_maxquant.py',
    '--distpath', dist_path,
    '--onedir',
    '--name', "test_2"])
