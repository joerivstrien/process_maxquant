"""
Author: Ariel Komen
Description: This code serves to construct a gui, select input files and run the process_maxquant script
"""
import sys
import logging
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QGroupBox, QFileDialog, \
    QVBoxLayout, QLineEdit, QHBoxLayout
from PyQt5.QtCore import pyqtSlot, QProcess

from process_maxquant import check_user_input
from process_maxquant import validate_user_parameters
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
        windowLayout = QHBoxLayout()
        windowLayout.addWidget(self.verticalGroupBox)
        windowLayout.addWidget(self.error_and_status_message_group_box)
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
        self.error_and_status_message_group_box = QGroupBox()
        main_vertical_layout = QVBoxLayout()
        error_and_status_vertical_layout = QVBoxLayout()

        maxquant_file_input_label = QLabel("Select the maxquant file: ")
        main_vertical_layout.addWidget(maxquant_file_input_label)

        status_label = QLabel("Program status: ")
        error_and_status_vertical_layout.addWidget(status_label)

        self.maxquant_file_input_field = QLineEdit()
        main_vertical_layout.addWidget(self.maxquant_file_input_field)

        self.status_message_label = QLabel("Program has not started")
        error_and_status_vertical_layout.addWidget(self.status_message_label)
        error_and_status_vertical_layout.addStretch(1)

        self.get_maxquant_file_button = QPushButton("Select maxquant file", self)
        self.get_maxquant_file_button.setToolTip("Select the maxquant file where the file should be in the csv format")
        self.get_maxquant_file_button.clicked.connect(self.openFileNameDialog)
        main_vertical_layout.addWidget(self.get_maxquant_file_button)

        settings_input_label = QLabel("Select the settings file: ")
        main_vertical_layout.addWidget(settings_input_label)

        error_status_label = QLabel("Error message: ")
        error_and_status_vertical_layout.addWidget(error_status_label)

        self.settings_file_input_field = QLineEdit()
        main_vertical_layout.addWidget(self.settings_file_input_field)

        self.error_message_label = QLabel("-")
        error_and_status_vertical_layout.addWidget(self.error_message_label)

        self.get_settings_file_button = QPushButton("Select the settings file", self)
        self.get_settings_file_button.setToolTip("Select the settings file where the file should be in the json format")
        self.get_settings_file_button.clicked.connect(self.openFileNameDialog)
        main_vertical_layout.addWidget(self.get_settings_file_button)

        self.process_maxquant_button = QPushButton('Process maxquant', self)
        self.process_maxquant_button.setToolTip('Clicking this button will start the program')
        self.process_maxquant_button.clicked.connect(self.execute_process_maxquant_script)
        main_vertical_layout.addWidget(self.process_maxquant_button)

        self.error_and_status_message_group_box.setLayout(error_and_status_vertical_layout)
        self.verticalGroupBox.setLayout(main_vertical_layout)

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

    def report_status(self, status_message):
        """"
        input:
        status_message = string
        output:
        None
        """
        self.status_message_label.setText(status_message)
        logging.info(status_message)
        QApplication.processEvents()  # This line causes the gui to display the message

    def report_error(self, error_message):
        """"
        Report the error message and close the process of the corresponding button.
        input:
        error_message = string
        output:
        None
        """
        self.error_message_label.setText(error_message)
        self.report_status("An error has occurred, please see the error report below.")

    @pyqtSlot()
    def execute_process_maxquant_script(self):
        self.process_maxquant_button.setEnabled(False)
        self.error_message_label.setText("-")
        logging.basicConfig(filename="process_maxquant_log.log", filemode="w", level=logging.DEBUG)
        if check_user_input(self, self.settings_file_input_field.text(), self.maxquant_file_input_field.text()) == False:
            self.process_maxquant_button.setEnabled(True)
            return

        settings_dict, is_json_loaded = load_json(self, self.settings_file_input_field.text())
        if is_json_loaded == False:
            self.process_maxquant_button.setEnabled(True)
            return

        protein_groups_dataframe, is_maxquant_file_loaded = read_in_protein_groups_file(self, self.maxquant_file_input_field.text())
        if is_maxquant_file_loaded == False:
            self.process_maxquant_button.setEnabled(True)
            return

        are_user_parameters_valid = validate_user_parameters(self, settings_dict, protein_groups_dataframe)
        if are_user_parameters_valid == False:
            self.process_maxquant_button.setEnabled(True)
            return

        protein_groups_dataframe, filtered_groups_dataframe, original_column_order = filter_dataframe_step(self, protein_groups_dataframe,
                                                                                    settings_dict)
        protein_groups_dataframe = fetch_uniprot_annotation_step(self, protein_groups_dataframe, settings_dict)

        protein_groups_dataframe = is_protein_in_mitocarta_step(self, settings_dict, protein_groups_dataframe)

        protein_groups_dataframe = apply_clustering_step(self, settings_dict, protein_groups_dataframe)

        dump_to_excel_step(self, protein_groups_dataframe, filtered_groups_dataframe, settings_dict, original_column_order)

        self.process_maxquant_button.setEnabled(True)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())