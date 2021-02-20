% Shortcut to open the main settings editor of goGPS
addPathGoGPS
log = Core.getLogger();
log.setOutMode(0,[],1);
ui = Core_UI.getInstance();
flag_wait = false;
ui.openGUI(flag_wait);
