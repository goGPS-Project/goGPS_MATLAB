% script to clean the workspace after a branch change
warning('off');
clear all
clear GO_Settings Logger Main_Settings Settings_Interface Command_Settings Core_UI Constellation_Collector Core Core_Sky ...
    Receiver_Work_Space Receiver_Output Receiver_Commons GNSS_Station Global_Configuration Command_Interpreter Core_Reference_Frame Meteo_Network Athmosphere Remote_Resource_Manager
% delete last_settings.ini
warning('off');
