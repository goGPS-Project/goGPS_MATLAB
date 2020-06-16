% script to clean the workspace after a branch change
warning('off');
GUI_Edit_Settings.closeGUI
clear all
clear Logger Prj_Settings Go_Settings Settings_Interface Command_Settings Core_UI GUI_Edit_Settings Core Remote_Resource_Manager Parallel_Manager Go_Slave Go_Wait_Bar Zernike
% delete last_settings.ini
warning('off');
