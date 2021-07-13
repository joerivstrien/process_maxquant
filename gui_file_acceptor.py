"""
Author: Ariel Komen
Description: This code serves to construct a gui, select input files and run the process_maxquant script
"""
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QGroupBox, QFileDialog, \
    QVBoxLayout, QLineEdit, QHBoxLayout
from PyQt5.QtCore import pyqtSlot

from process_maxquant import load_json
from process_maxquant import read_in_protein_groups_file
from process_maxquant import filter_dataframe_step
from process_maxquant import fetch_uniprot_annotation_step
from process_maxquant import is_protein_in_mitocarta_step
from process_maxquant import apply_clustering_step
from process_maxquant import dump_to_excel_step
#ToDo, Make sure that the error handling is applied in a way that makes sense for gui's.

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
        """"
        Create a gui and put each component in a vertical layout
        input:
        None
        output:
        None
        """
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.verticalGroupBox = QGroupBox()
        vertical_layout = QVBoxLayout()
        horizontal_status_layout = QHBoxLayout()
        horizontal_message_layout = QHBoxLayout()
        horizontal_error_status_layout = QHBoxLayout()
        horizontal_error_message_layout = QHBoxLayout()

        maxquant_file_input_label = QLabel("Select the maxquant file: ")
        status_label = QLabel("Program status: ")
        horizontal_status_layout.addWidget(maxquant_file_input_label)
        horizontal_status_layout.addStretch(1)
        horizontal_status_layout.addWidget(status_label, 1)
        vertical_layout.addLayout(horizontal_status_layout)

        self.maxquant_file_input_field = QLineEdit()
        self.status_message_label = QLabel("Program has not started")
        horizontal_message_layout.addWidget(self.maxquant_file_input_field)
        horizontal_message_layout.addStretch(1)
        horizontal_message_layout.addWidget(self.status_message_label, 1)
        vertical_layout.addLayout(horizontal_message_layout)

        self.get_maxquant_file_button = QPushButton("Select maxquant file", self)
        self.get_maxquant_file_button.setToolTip("Select the maxquant file where the file should be in the csv format")
        self.get_maxquant_file_button.clicked.connect(self.openFileNameDialog)
        vertical_layout.addWidget(self.get_maxquant_file_button)

        settings_input_label = QLabel("Select the settings file: ")
        error_status_label = QLabel("Error message: ")
        horizontal_error_status_layout.addWidget(settings_input_label)
        horizontal_error_status_layout.addStretch(1)
        horizontal_error_status_layout.addWidget(error_status_label, 1)
        vertical_layout.addLayout(horizontal_error_status_layout)

        self.settings_file_input_field = QLineEdit()
        self.error_message_label = QLabel("-")
        horizontal_error_message_layout.addWidget(self.settings_file_input_field)
        horizontal_error_message_layout.addStretch(1)
        horizontal_error_message_layout.addWidget(self.error_message_label, 1)
        vertical_layout.addLayout(horizontal_error_message_layout)

        self.get_settings_file_button = QPushButton("Select the settings file", self)
        self.get_settings_file_button.setToolTip("Select the settings file where the file should be in the json format")
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

    def report_status(self):
        pass

    def report_error(self):
        pass

    @pyqtSlot()
    def execute_process_maxquant_script(self):
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