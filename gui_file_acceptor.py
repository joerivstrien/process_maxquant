"""
Author: Ariel Komen
Description: This code serves to make a pyqt5 gui
"""
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QGroupBox, QFileDialog, \
    QVBoxLayout, QLineEdit
from PyQt5.QtCore import pyqtSlot

from process_maxquant import load_json
from process_maxquant import read_in_protein_groups_file
from process_maxquant import filter_dataframe_step
from process_maxquant import fetch_uniprot_annotation_step
from process_maxquant import is_protein_in_mitocarta_step
from process_maxquant import apply_clustering_step
from process_maxquant import dump_to_excel_step

class App(QWidget):
    def __init__(self):
        super().__init__()
        self.title = 'Process maxquant program'
        self.left = 10
        self.top = 10
        self.width = 500
        self.height = 400

        self.create_gui_layout()
        windowLayout = QVBoxLayout()
        windowLayout.addWidget(self.verticalGroupBox)
        self.setLayout(windowLayout)

        self.show()

    def create_gui_layout(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.verticalGroupBox = QGroupBox()
        vertical_layout = QVBoxLayout()

        maxquant_file_input_label = QLabel("Select the maxquant file: ")
        vertical_layout.addWidget(maxquant_file_input_label)

        self.maxquant_file_input_field = QLineEdit()
        vertical_layout.addWidget(self.maxquant_file_input_field)

        self.get_maxquant_file_button = QPushButton("Select maxquant file", self)
        self.get_maxquant_file_button.setToolTip("Select the maxquant file in csv format")
        self.get_maxquant_file_button.clicked.connect(self.openFileNameDialog)
        vertical_layout.addWidget(self.get_maxquant_file_button)

        settings_input_label = QLabel("Select the settings file: ")
        vertical_layout.addWidget(settings_input_label)

        self.settings_file_input_field = QLineEdit()
        vertical_layout.addWidget(self.settings_file_input_field)

        self.get_settings_file_button = QPushButton("Select the settings file", self)
        self.get_settings_file_button.setToolTip("Select the settings file in the json format")
        self.get_settings_file_button.clicked.connect(self.openFileNameDialog)
        vertical_layout.addWidget(self.get_settings_file_button)

        process_maxquant_button = QPushButton('Process maxquant', self)
        process_maxquant_button.setToolTip('Clicking this button will start the program')
        process_maxquant_button.clicked.connect(self.execute_process_maxquant_script)
        vertical_layout.addWidget(process_maxquant_button)

        self.verticalGroupBox.setLayout(vertical_layout)

    @pyqtSlot()
    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Python Files (*.py)", options=options)
        if self.sender() is self.get_settings_file_button:
            self.settings_file_input_field.setText(file_name)
        elif self.sender() is self.get_maxquant_file_button:
            self.maxquant_file_input_field.setText(file_name)

    @pyqtSlot()
    def execute_process_maxquant_script(self):
        #do all the stuff I would do in process_maxquant.
        settings_dict = load_json(self.settings_file_input_field.text())
        protein_groups_dataframe = read_in_protein_groups_file(self.maxquant_file_input_field.text())

        protein_groups_dataframe, filtered_groups_dataframe = filter_dataframe_step(protein_groups_dataframe,
                                                                                    settings_dict)
        protein_groups_dataframe = fetch_uniprot_annotation_step(protein_groups_dataframe, settings_dict)

        protein_groups_dataframe = is_protein_in_mitocarta_step(settings_dict, protein_groups_dataframe)

        protein_groups_dataframe = apply_clustering_step(settings_dict, protein_groups_dataframe)

        dump_to_excel_step(protein_groups_dataframe, filtered_groups_dataframe, settings_dict)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())