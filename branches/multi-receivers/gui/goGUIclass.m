classdef goGUIclass < handle

    properties (Constant)
        % Constellations id => naming inside the goObservation object
        isWin = 0;
        isUnix = 1;
        %isLinux = 2;
        %isMac = 3;
        isRinex = 0;
        isBin = 1;
    end

    
    properties (GetAccess = 'private', SetAccess = 'private')
    
    % =========================================================================
    %   Handlers
    % =========================================================================
    
        goh = [];                 % goGPS gui handler
        
        interfaceOS = 0;          % 0 = Windows
                                  % 1 = Linux
                                  % 2 = Mac
                                  
        % structure containing the state of the parameters of the figure
        status; 
    end
    
    % =========================================================================
    %    PUBLIC METHODS
    % =========================================================================

    methods
        % Creator (Brahma)
        function obj = goGUIclass(handles, typeOS)
            if (length(intersect(typeOS, [obj.isWin obj.isUnix])) ~= 1)
                typeOS = obj.isWin;
            end
            obj.interfaceOS = typeOS;
            
            obj.init(handles);
        end
        
        % Destructor (Shiva)
        function delete(obj)
        end

    % =========================================================================
    %    GETTER
    % =========================================================================

        % Get the mode as integer value
        function mode = getMode(obj)
            contents_mode = cellstr(get(obj.goh.mode,'String'));
            contents_nav_mon = cellstr(get(obj.goh.nav_mon,'String'));
            contents_kalman_ls = cellstr(get(obj.goh.kalman_ls,'String'));
            contents_code_dd_sa = cellstr(get(obj.goh.code_dd_sa,'String'));
            if (strcmp(contents_mode{get(obj.goh.mode,'Value')},'Post-processing'))
                if (strcmp(contents_kalman_ls{get(obj.goh.kalman_ls,'Value')},'Kalman filter'))
                    if (strcmp(contents_code_dd_sa{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference'))
                        mode = 14;
                    elseif (strcmp(contents_code_dd_sa{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference for several receivers'))
                        mode = 15;
                    elseif (strcmp(contents_code_dd_sa{get(obj.goh.code_dd_sa,'Value')},'Code and phase stand-alone'))
                        mode = 4;
                    elseif (strcmp(contents_code_dd_sa{get(obj.goh.code_dd_sa,'Value')},'Code double difference'))
                        mode = 12;
                    else %Code stand-alone
                        mode = 2;
                    end
                else %Least squares
                    if (strcmp(contents_code_dd_sa{get(obj.goh.code_dd_sa,'Value')},'Code double difference'))
                        mode = 11;
                    else %Code stand-alone
                        mode = 1;
                    end
                end
            else %Real-time
                if (strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Navigation'))
                    mode = 24;
                elseif (strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Rover monitor'))
                    mode = 21;
                elseif (strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Master monitor'))
                    mode = 22;
                else %Rover and Master monitor
                    mode = 23;
                end
            end
        end
        
    % =========================================================================
    %    OUTPUT FUNCTION
    % =========================================================================
    
        function funout = outputFun(obj)        
            mode = obj.getMode();
            mode_vinc = get(obj.goh.constraint,'Value');
            if (get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files)
                mode_data = 0;
            else %goGPS data
                mode_data = 1;
            end
            contents_dyn_mod = cellstr(get(obj.goh.dyn_mod,'String'));
            if (strcmp(contents_dyn_mod{get(obj.goh.dyn_mod,'Value')},'Variable') || get(obj.goh.stopGOstop,'Value'))
                flag_var_dyn_model = 1;
            else
                flag_var_dyn_model = 0;
            end
            mode_ref = get(obj.goh.ref_path,'Value');
            flag_ms_pos = get(obj.goh.master_pos,'Value');
            flag_ms = get(obj.goh.plot_master,'Value');
            flag_ge = get(obj.goh.google_earth,'Value');
            flag_cov = get(obj.goh.err_ellipse,'Value');
            flag_NTRIP = get(obj.goh.use_ntrip,'Value');
            flag_amb = get(obj.goh.plot_amb,'Value');
            flag_skyplot = get(obj.goh.no_skyplot_snr,'Value');
            flag_plotproc = get(obj.goh.plotproc,'Value');
            flag_stopGOstop = get(obj.goh.stopGOstop,'Value');
            filerootIN = get(obj.goh.gogps_data_input,'String');
            filerootOUT = [get(obj.goh.gogps_data_output,'String') '\' get(obj.goh.gogps_data_output_prefix,'String')];
            filerootIN(filerootIN == '\') = '/';
            filerootOUT(filerootOUT == '\') = '/';
            i = 1;
            j = length(filerootOUT);
            while (~isempty(dir([filerootOUT '_????_rover.bin'])) || ...
                    ~isempty(dir([filerootOUT '_master*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_????_obs*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_????_eph*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_????_dyn*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_sat*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_kal*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_dt*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_conf*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_dop*.bin'])) || ...
                    ~isempty(dir([filerootOUT '_ECEF*.txt'])) || ...
                    ~isempty(dir([filerootOUT '_geod*.txt'])) || ...
                    ~isempty(dir([filerootOUT '_plan*.txt'])) || ...
                    ~isempty(dir([filerootOUT '_????_NMEA*.txt'])) || ...
                    ~isempty(dir([filerootOUT '.kml'])) )
                
                filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
                i = i + 1;
            end
            filename_R_obs = get(obj.goh.RINEX_rover_obs,'String');
            filename_M_obs = get(obj.goh.RINEX_master_obs,'String');
            filename_nav = get(obj.goh.RINEX_nav,'String');
            if (strcmpi(filename_nav(end-3:end),'.sp3'))
                flag_SP3 = 1;
            else
                flag_SP3 = 0;
            end
            filename_ref = get(obj.goh.ref_path_input,'String');
            
            contents = cellstr(get(obj.goh.crs,'String'));
            if (strcmp(contents{get(obj.goh.crs,'Value')},'ECEF (X,Y,Z)'))
                XM = str2double(get(obj.goh.master_X,'String'));
                YM = str2double(get(obj.goh.master_Y,'String'));
                ZM = str2double(get(obj.goh.master_Z,'String'));
            else
                latM = str2double(get(obj.goh.master_lat,'String'));
                lonM = str2double(get(obj.goh.master_lon,'String'));
                hM = str2double(get(obj.goh.master_h,'String'));
                [XM, YM, ZM] = geod2cart (latM*pi/180, lonM*pi/180, hM, 6378137, 1/298.257222101);
            end
            pos_M_man = [XM; YM; ZM];
            
            contents = cellstr(get(obj.goh.num_receivers,'String'));
            num_rec = str2double(contents{get(obj.goh.num_receivers,'Value')});
            
            if num_rec >= 1
                contentsProt = cellstr(get(obj.goh.protocol_select_0,'String'));
                if (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'UBX (u-blox)'))
                    protocol_idx(1) = 0;
                elseif (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'iTalk (Fastrax)'))
                    protocol_idx(1) = 1;
                elseif (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'SkyTraq'))
                    protocol_idx(1) = 2;
                end
                
                if num_rec >= 2
                    contentsProt = cellstr(get(obj.goh.protocol_select_1,'String'));
                    if (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'UBX (u-blox)'))
                        protocol_idx(2) = 0;
                    elseif (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'iTalk (Fastrax)'))
                        protocol_idx(2) = 1;
                    elseif (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'SkyTraq'))
                        protocol_idx(2) = 2;
                    end
                    
                    if num_rec >= 3
                        contentsProt = cellstr(get(obj.goh.protocol_select_2,'String'));
                        if (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'UBX (u-blox)'))
                            protocol_idx(3) = 0;
                        elseif (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'iTalk (Fastrax)'))
                            protocol_idx(3) = 1;
                        elseif (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'SkyTraq'))
                            protocol_idx(3) = 2;
                        end
                        
                        if num_rec >= 4
                            contentsProt = cellstr(get(obj.goh.protocol_select_3,'String'));
                            if (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'UBX (u-blox)'))
                                protocol_idx(4) = 0;
                            elseif (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'iTalk (Fastrax)'))
                                protocol_idx(4) = 1;
                            elseif (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'SkyTraq'))
                                protocol_idx(4) = 2;
                            end
                        end
                    end
                end
            end
            
            funout = cell(23,1);
            
            funout{1} = mode;
            funout{2} = mode_vinc;
            funout{3} = mode_data;
            funout{4} = mode_ref;
            funout{5} = flag_ms_pos;
            funout{6} = flag_ms;
            funout{7} = flag_ge;
            funout{8} = flag_cov;
            funout{9} = flag_NTRIP;
            funout{10} = flag_amb;
            funout{11} = flag_skyplot;
            funout{12} = flag_plotproc;
            funout{13} = flag_var_dyn_model;
            funout{14} = flag_stopGOstop;
            funout{15} = flag_SP3;
            funout{16} = filerootIN;
            funout{17} = filerootOUT;
            funout{18} = filename_R_obs;
            funout{19} = filename_M_obs;
            funout{20} = filename_nav;
            funout{21} = filename_ref;
            funout{22} = pos_M_man;
            funout{23} = protocol_idx;
            
            global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
            global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
            global min_nsat cutoff snr_threshold cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
            global h_antenna
            global tile_header tile_georef dtm_dir
            global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
            global nmea_init
            global flag_doppler_cs
            global COMportR
            
            contents = cellstr(get(obj.goh.com_select_0,'String'));
            COMportR0 = contents{get(obj.goh.com_select_0,'Value')};
            contents = cellstr(get(obj.goh.com_select_1,'String'));
            COMportR1 = contents{get(obj.goh.com_select_1,'Value')};
            contents = cellstr(get(obj.goh.com_select_2,'String'));
            COMportR2 = contents{get(obj.goh.com_select_2,'Value')};
            contents = cellstr(get(obj.goh.com_select_3,'String'));
            COMportR3 = contents{get(obj.goh.com_select_3,'Value')};
            
            if num_rec >= 1
                COMportR{1,1} = COMportR0;
                if num_rec >= 2
                    COMportR{2,1} = COMportR1;
                    if num_rec >= 3
                        COMportR{3,1} = COMportR2;
                        if num_rec >= 4
                            COMportR{4,1} = COMportR3;
                        end
                    end
                end
            end
            
            flag_doppler_cs = get(obj.goh.flag_doppler,'Value');
            sigmaq0 = str2double(get(obj.goh.std_init,'String'))^2;
            sigmaq_vE = str2double(get(obj.goh.std_X,'String'))^2;
            sigmaq_vN = str2double(get(obj.goh.std_Y,'String'))^2;
            sigmaq_vU = str2double(get(obj.goh.std_Z,'String'))^2;
            sigmaq_vel = str2double(get(obj.goh.std_vel,'String'))^2;
            sigmaq_cod1 = str2double(get(obj.goh.std_code,'String'))^2;
            sigmaq_cod2 = 0.16;
            if (get(obj.goh.toggle_std_phase,'Value'))
                sigmaq_ph = str2double(get(obj.goh.std_phase,'String'))^2;
            else
                sigmaq_ph = 1e30;
            end
            sigmaq0_N = 100;
            if (get(obj.goh.toggle_std_dtm,'Value'))
                sigmaq_dtm = str2double(get(obj.goh.std_dtm,'String'))^2;
            else
                sigmaq_dtm = 1e30;
            end
            min_nsat = str2double(get(obj.goh.min_sat,'String'));
            if (mode == 2)
                disp('Minimum number of satellites is forced to 4 (for stand-alone positioning)');
                min_nsat = 4;
            end
            cutoff = str2double(get(obj.goh.cut_off,'String'));
            snr_threshold = str2double(get(obj.goh.snr_thres,'String'));
            cs_threshold = str2double(get(obj.goh.cs_thresh,'String'));
            if (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_0)
                weights = 0;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_1)
                weights = 1;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_2)
                weights = 2;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_3)
                weights = 3;
            end
            snr_a = 30;
            snr_0 = 10;
            snr_1 = 50;
            snr_A = 30;
            obj.selectAmbiguityRestartMethod();
            obj.selectDynMode();
            o1 = order;
            o2 = order*2;
            o3 = order*3;
            h_antenna = str2double(get(obj.goh.antenna_h,'String'));
            dtm_dir = get(obj.goh.dtm_path,'String');
            try
                load([dtm_dir '/tiles/tile_header'], 'tile_header');
                load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
            catch
                tile_header.nrows = 0;
                tile_header.ncols = 0;
                tile_header.cellsize = 0;
                tile_header.nodata = 0;
                tile_georef = zeros(1,1,4);
            end
            master_ip = get(obj.goh.IP_address,'String');
            master_port = str2double(get(obj.goh.port,'String'));
            ntrip_user = get(obj.goh.username,'String');
            ntrip_pw = get(obj.goh.password,'Userdata');
            ntrip_mountpoint = get(obj.goh.mountpoint,'String');
            phiApp = str2double(get(obj.goh.approx_lat,'String'));
            lamApp = str2double(get(obj.goh.approx_lon,'String'));
            hApp = str2double(get(obj.goh.approx_h,'String'));
            [XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, 6378137, 1/298.257222101);
            if ~isnan(XApp) && ~isnan(YApp) && ~isnan(ZApp)
                nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);
            else
                nmea_init = '';
            end
        end
    % =========================================================================
    %    INTERFACE MANAGEMENT
    % =========================================================================

        %   PASSWORD
        % =================================================================
        
        % Show the password in the password field
        function showPassword(obj)
            if get(obj.goh.show_password,'Value')
                set(obj.goh.password,'String',obj.status.password);
            else
                SizePass = size(obj.status.password); % Find the number of asterisks
                if SizePass(2) > 0
                    asterisk(1,1:SizePass(2)) = '*'; % Create a string of asterisks the same size as the password
                    set(obj.goh.password,'String',asterisk) % Set the text in the password edit box to the asterisk string
                else
                    set(obj.goh.password,'String','')
                end

            end
        end
                
        % If the showPassword flag is enabled set a visible password
        % otherwise save the new password and show asterisks
        function setPassword(obj, password)
            if (nargin == 2)
                obj.status.password = password;
            end
            obj.showPassword()
        end
        
        % Get the password
        function password = getPassword(obj)
            password = obj.status.password;
        end
        
        % Modify the password
        function modifyPassword(obj,newkey, newchar)
            password = obj.getPassword();
            switch newkey
                case 'backspace'
                    password = password(1:end-1); % Delete the last character in the password
                case 'delete'
                    password = password(2:end); % Delete the first character in the password
                otherwise
                    % If pressed key produces a printable character
                    if (uint8(newchar) > 32)
                        password = [password newchar]; % Add the typed character to the password
                        pause(0.001)%to avoid unwanted character output before the cursor
                    end
            end
            
            obj.setPassword(password); % Store the password in its current state
        end
        
        %   TYPE SELECTORS
        % =================================================================

        % Select the kind of approach for the computation (real time - post proc)
        function selectProcassingApproach(obj)
            contents = cellstr(get(obj.goh.mode,'String'));
            if (strcmp(contents{get(obj.goh.mode,'Value')},'Real-time'))
                try
                    instrhwinfo;
                catch
                    warndlg('Instrument Control Toolbox is needed to run goGPS in real-time mode.', 'Warning');
                    set(obj.goh.mode, 'Value', 2);
                    obj.selectProcassingApproach();
                    return
                end
                set(obj.goh.nav_mon, 'Enable', 'on');
                set(obj.goh.kalman_ls, 'Enable', 'off');
                set(obj.goh.kalman_ls, 'Value', 1);
                set(obj.goh.code_dd_sa, 'Enable', 'off');
                set(obj.goh.code_dd_sa, 'Value', 4);
                set(obj.goh.rinex_files, 'Enable', 'off');
                set(obj.goh.gogps_data, 'Enable', 'off');
                
                set(obj.goh.plot_amb, 'Enable', 'off');
                set(obj.goh.no_skyplot_snr, 'Enable', 'on');
                set(obj.goh.plotproc, 'Enable', 'on');
                set(obj.goh.flag_doppler, 'Enable', 'off');
                obj.togglePlotOnProcessing();
                
                obj.selectReceiverMode();
                
                %disable file input fields
                set(obj.goh.RINEX_rover_obs, 'Enable', 'off');
                set(obj.goh.RINEX_master_obs, 'Enable', 'off');
                set(obj.goh.RINEX_nav, 'Enable', 'off');
                set(obj.goh.browse_rover_obs, 'Enable', 'off');
                set(obj.goh.browse_master_obs, 'Enable', 'off');
                set(obj.goh.browse_nav, 'Enable', 'off');
                set(obj.goh.text_RINEX_rover_obs, 'Enable', 'off');
                set(obj.goh.text_RINEX_master_obs, 'Enable', 'off');
                set(obj.goh.text_RINEX_nav, 'Enable', 'off');
                set(obj.goh.gogps_data_input, 'Enable', 'off');
                set(obj.goh.browse_gogps_input, 'Enable', 'off');
                set(obj.goh.text_gogps_input, 'Enable', 'off');
            else
                set(obj.goh.nav_mon, 'Enable', 'off');
                set(obj.goh.nav_mon, 'Value', 1);
                
                obj.selectReceiverMode();

                set(obj.goh.kalman_ls, 'Enable', 'on');
                set(obj.goh.code_dd_sa, 'Enable', 'on');
                set(obj.goh.rinex_files, 'Enable', 'on');
                set(obj.goh.gogps_data, 'Enable', 'on');
                set(obj.goh.text_num_receivers, 'Enable', 'off');
                set(obj.goh.num_receivers, 'Enable', 'off');
                set(obj.goh.com_select_0, 'Enable', 'off');
                set(obj.goh.com_select_1, 'Enable', 'off');
                set(obj.goh.com_select_2, 'Enable', 'off');
                set(obj.goh.com_select_3, 'Enable', 'off');
                set(obj.goh.protocol_select_0, 'Enable', 'off');
                set(obj.goh.protocol_select_1, 'Enable', 'off');
                set(obj.goh.protocol_select_2, 'Enable', 'off');
                set(obj.goh.protocol_select_3, 'Enable', 'off');
                set(obj.goh.use_ntrip, 'Enable', 'off');
                contents = cellstr(get(obj.goh.code_dd_sa,'String'));
                if (get(obj.goh.plotproc,'Value') && (strcmp(contents{get(obj.goh.code_dd_sa,'Value')}, ...
                        'Code and phase double difference') || strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase stand-alone') || ...
                        strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference for several receivers')))
                    set(obj.goh.plot_amb, 'Enable', 'on');
                    obj.toggleAmbiguities();
                end
                
                %enable/disable file input fields
                if(get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files);
                    obj.selectFileType(obj.goh.rinex_files);
                else
                    obj.selectFileType(obj.goh.gogps_data);
                end
                
                obj.selectProcessingType();
                
                %disable approximate position
                set(obj.goh.text_approx_pos, 'Enable', 'off');
                set(obj.goh.approx_lat, 'Enable', 'off');
                set(obj.goh.approx_lon, 'Enable', 'off');
                set(obj.goh.approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat, 'Enable', 'off');
                set(obj.goh.text_approx_lon, 'Enable', 'off');
                set(obj.goh.text_approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat_unit, 'Enable', 'off');
                set(obj.goh.text_approx_lon_unit, 'Enable', 'off');
                set(obj.goh.text_approx_h_unit, 'Enable', 'off');
                
                %disable NTRIP parameters
                set(obj.goh.IP_address, 'Enable', 'off');
                set(obj.goh.port, 'Enable', 'off');
                set(obj.goh.mountpoint, 'Enable', 'off');
                set(obj.goh.username, 'Enable', 'off');
                set(obj.goh.password, 'Enable', 'off');
                set(obj.goh.show_password, 'Enable', 'off');
                set(obj.goh.text_IP_address, 'Enable', 'off');
                set(obj.goh.text_port, 'Enable', 'off');
                set(obj.goh.text_mountpoint, 'Enable', 'off');
                set(obj.goh.text_username, 'Enable', 'off');
                set(obj.goh.text_password, 'Enable', 'off');
                
                obj.toggleStopGoStop();
            end            
        end
        
        % Select navigation / monitor mode
        function selectReceiverMode(obj)
            contents = cellstr(get(obj.goh.nav_mon,'String'));
            if (strcmp(contents{get(obj.goh.nav_mon,'Value')},'Navigation'))
                
                %enable options
                set(obj.goh.constraint, 'Enable', 'on');
                set(obj.goh.ref_path, 'Enable', 'on');
                obj.toggleReferencePath();
                set(obj.goh.plot_master, 'Enable', 'on');
                set(obj.goh.err_ellipse, 'Enable', 'on');
                set(obj.goh.google_earth, 'Enable', 'on');
                set(obj.goh.text_num_receivers, 'Enable', 'off');
                set(obj.goh.num_receivers, 'Value', 1);
                set(obj.goh.num_receivers, 'Enable', 'off');
                set(obj.goh.com_select_0, 'Enable', 'on');
                set(obj.goh.com_select_1, 'Enable', 'off');
                set(obj.goh.com_select_2, 'Enable', 'off');
                set(obj.goh.com_select_3, 'Enable', 'off');
                set(obj.goh.protocol_select_0, 'Enable', 'on');
                set(obj.goh.protocol_select_1, 'Enable', 'off');
                set(obj.goh.protocol_select_2, 'Enable', 'off');
                set(obj.goh.protocol_select_3, 'Enable', 'off');
                set(obj.goh.use_ntrip, 'Enable', 'on');
                set(obj.goh.no_skyplot_snr, 'Enable', 'on');
                
                set(obj.goh.gogps_data_output, 'Enable', 'on');
                set(obj.goh.text_gogps_data_output, 'Enable', 'on');
                set(obj.goh.browse_gogps_data_output, 'Enable', 'on');
                set(obj.goh.gogps_data_output_prefix, 'Enable', 'on');
                set(obj.goh.text_gogps_data_output_prefix, 'Enable', 'on');
                
                obj.selectProcessingType();
                
                set(obj.goh.plotproc, 'Enable', 'on');
                obj.togglePlotOnProcessing();
                
                set(obj.goh.stopGOstop, 'Enable', 'on');
                set(obj.goh.text_stopGOstop, 'Enable', 'on');
                obj.toggleStopGoStop();
                
                cell_contents = cell(4,1);
                cell_contents{1} = 'Const. velocity';
                cell_contents{2} = 'Const. acceleration';
                cell_contents{3} = 'Static';
                cell_contents{4} = 'Variable';
                set(obj.goh.dyn_mod, 'String', cell_contents);
                
                %enable weights
                set(obj.goh.weight_0, 'Enable', 'on');
                set(obj.goh.weight_1, 'Enable', 'on');
                set(obj.goh.weight_2, 'Enable', 'on');
                set(obj.goh.weight_3, 'Enable', 'on');
                
                %enable master connection
                set(obj.goh.IP_address, 'Enable', 'on');
                set(obj.goh.port, 'Enable', 'on');
                set(obj.goh.text_IP_address, 'Enable', 'on');
                set(obj.goh.text_port, 'Enable', 'on');
                obj.toggleNTRIP();                
                %disable approximate position
                set(obj.goh.text_approx_pos, 'Enable', 'off');
                set(obj.goh.approx_lat, 'Enable', 'off');
                set(obj.goh.approx_lon, 'Enable', 'off');
                set(obj.goh.approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat, 'Enable', 'off');
                set(obj.goh.text_approx_lon, 'Enable', 'off');
                set(obj.goh.text_approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat_unit, 'Enable', 'off');
                set(obj.goh.text_approx_lon_unit, 'Enable', 'off');
                set(obj.goh.text_approx_h_unit, 'Enable', 'off');
                
            else
                
                %disable some options
                set(obj.goh.constraint, 'Enable', 'off');
                set(obj.goh.ref_path, 'Enable', 'off');
                set(obj.goh.plot_master, 'Enable', 'off');
                set(obj.goh.err_ellipse, 'Enable', 'off');
                set(obj.goh.google_earth, 'Enable', 'off');
                set(obj.goh.no_skyplot_snr, 'Enable', 'off');
                set(obj.goh.plotproc, 'Enable', 'off');
                set(obj.goh.plot_amb, 'Enable', 'off');
                %set(obj.goh.gogps_data_output, 'Enable', 'off');
                %set(obj.goh.text_gogps_data_output, 'Enable', 'off');
                %set(obj.goh.browse_gogps_data_output, 'Enable', 'off');
                %set(obj.goh.gogps_data_output_prefix, 'Enable', 'off');
                %set(obj.goh.text_gogps_data_output_prefix, 'Enable', 'off');
                set(obj.goh.ref_path_input, 'Enable', 'off');
                set(obj.goh.text_ref_path_input, 'Enable', 'off');
                set(obj.goh.browse_ref_path_input, 'Enable', 'off');
                
                %disable Kalman filters settings
                set(obj.goh.std_X, 'Enable', 'off');
                set(obj.goh.std_Y, 'Enable', 'off');
                set(obj.goh.std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X, 'Enable', 'off');
                set(obj.goh.text_std_Y, 'Enable', 'off');
                set(obj.goh.text_std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X_unit, 'Enable', 'off');
                set(obj.goh.text_std_Y_unit, 'Enable', 'off');
                set(obj.goh.text_std_Z_unit, 'Enable', 'off');
                set(obj.goh.std_code, 'Enable', 'off');
                set(obj.goh.std_phase, 'Enable', 'off');
                set(obj.goh.text_std_code, 'Enable', 'off');
                set(obj.goh.toggle_std_phase, 'Enable', 'off');
                set(obj.goh.text_std_code_unit, 'Enable', 'off');
                set(obj.goh.text_std_phase_unit, 'Enable', 'off');
                set(obj.goh.std_init, 'Enable', 'off');
                set(obj.goh.std_dtm, 'Enable', 'off');
                set(obj.goh.std_vel, 'Enable', 'off');
                set(obj.goh.text_std_init, 'Enable', 'off');
                set(obj.goh.toggle_std_dtm, 'Enable', 'off');
                set(obj.goh.text_std_vel, 'Enable', 'off');
                set(obj.goh.text_std_init_unit, 'Enable', 'off');
                set(obj.goh.text_std_dtm_unit, 'Enable', 'off');
                set(obj.goh.text_std_vel_unit, 'Enable', 'off');
                set(obj.goh.cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh_unit, 'Enable', 'off');
                set(obj.goh.cut_off, 'Enable', 'off');
                set(obj.goh.text_cut_off, 'Enable', 'off');
                set(obj.goh.text_cut_off_unit, 'Enable', 'off');
                set(obj.goh.snr_thres, 'Enable', 'off');
                set(obj.goh.text_snr_thres, 'Enable', 'off');
                set(obj.goh.text_snr_thres_unit, 'Enable', 'off');
                set(obj.goh.min_sat, 'Enable', 'off');
                set(obj.goh.text_min_sat, 'Enable', 'off');
                set(obj.goh.antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h_unit, 'Enable', 'off');
                set(obj.goh.dtm_path, 'Enable', 'off');
                set(obj.goh.text_dtm_path, 'Enable', 'off');
                set(obj.goh.browse_dtm_path, 'Enable', 'off');
                set(obj.goh.weight_0, 'Enable', 'off');
                set(obj.goh.weight_1, 'Enable', 'off');
                set(obj.goh.weight_2, 'Enable', 'off');
                set(obj.goh.weight_3, 'Enable', 'off');
                set(obj.goh.master_pos, 'Enable', 'off');
                set(obj.goh.crs, 'Enable', 'off');
                set(obj.goh.master_X, 'Enable', 'off');
                set(obj.goh.master_Y, 'Enable', 'off');
                set(obj.goh.master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X, 'Enable', 'off');
                set(obj.goh.text_master_Y, 'Enable', 'off');
                set(obj.goh.text_master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X_unit, 'Enable', 'off');
                set(obj.goh.text_master_Y_unit, 'Enable', 'off');
                set(obj.goh.text_master_Z_unit, 'Enable', 'off');
                set(obj.goh.master_lat, 'Enable', 'off');
                set(obj.goh.master_lon, 'Enable', 'off');
                set(obj.goh.master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat, 'Enable', 'off');
                set(obj.goh.text_master_lon, 'Enable', 'off');
                set(obj.goh.text_master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat_unit, 'Enable', 'off');
                set(obj.goh.text_master_lon_unit, 'Enable', 'off');
                set(obj.goh.text_master_h_unit, 'Enable', 'off');
                set(obj.goh.flag_doppler, 'Enable', 'off');
                set(obj.goh.amb_select, 'Enable', 'off');
                
                cell_contents = cell(2,1);
                cell_contents{1} = 'Constant';
                cell_contents{2} = 'Variable';
                old_value = get(obj.goh.dyn_mod, 'Value');
                if (old_value == 3), set(obj.goh.dyn_mod, 'Value', 1); end
                if (old_value == 4), set(obj.goh.dyn_mod, 'Value', 2); end
                set(obj.goh.dyn_mod, 'String', cell_contents);
                
                if (strcmp(contents{get(obj.goh.nav_mon,'Value')},'Rover monitor'))
                    
                    obj.setNumReceiverIn();
                    set(obj.goh.use_ntrip, 'Enable', 'off');
                    
                    %disable approximate position
                    set(obj.goh.text_approx_pos, 'Enable', 'off');
                    set(obj.goh.approx_lat, 'Enable', 'off');
                    set(obj.goh.approx_lon, 'Enable', 'off');
                    set(obj.goh.approx_h, 'Enable', 'off');
                    set(obj.goh.text_approx_lat, 'Enable', 'off');
                    set(obj.goh.text_approx_lon, 'Enable', 'off');
                    set(obj.goh.text_approx_h, 'Enable', 'off');
                    set(obj.goh.text_approx_lat_unit, 'Enable', 'off');
                    set(obj.goh.text_approx_lon_unit, 'Enable', 'off');
                    set(obj.goh.text_approx_h_unit, 'Enable', 'off');
                    
                    %disable NTRIP parameters
                    set(obj.goh.IP_address, 'Enable', 'off');
                    set(obj.goh.port, 'Enable', 'off');
                    set(obj.goh.mountpoint, 'Enable', 'off');
                    set(obj.goh.username, 'Enable', 'off');
                    set(obj.goh.password, 'Enable', 'off');
                    set(obj.goh.show_password, 'Enable', 'off');
                    set(obj.goh.text_IP_address, 'Enable', 'off');
                    set(obj.goh.text_port, 'Enable', 'off');
                    set(obj.goh.text_mountpoint, 'Enable', 'off');
                    set(obj.goh.text_username, 'Enable', 'off');
                    set(obj.goh.text_password, 'Enable', 'off');
                    
                    set(obj.goh.stopGOstop, 'Enable', 'on');
                    set(obj.goh.text_stopGOstop, 'Enable', 'on');
                    obj.toggleStopGoStop();
                    
                elseif (strcmp(contents{get(obj.goh.nav_mon,'Value')},'Master monitor'))
                    
                    set(obj.goh.text_num_receivers, 'Enable', 'off');
                    set(obj.goh.num_receivers, 'Enable', 'off');
                    set(obj.goh.com_select_0, 'Enable', 'off');
                    set(obj.goh.com_select_1, 'Enable', 'off');
                    set(obj.goh.com_select_2, 'Enable', 'off');
                    set(obj.goh.com_select_3, 'Enable', 'off');
                    set(obj.goh.protocol_select_0, 'Enable', 'off');
                    set(obj.goh.protocol_select_1, 'Enable', 'off');
                    set(obj.goh.protocol_select_2, 'Enable', 'off');
                    set(obj.goh.protocol_select_3, 'Enable', 'off');
                    set(obj.goh.use_ntrip, 'Enable', 'on');
                    set(obj.goh.dyn_mod, 'Enable', 'off');
                    %set(obj.goh.text_dyn_mod, 'Enable', 'off');
                    
                    %enable master connection
                    set(obj.goh.IP_address, 'Enable', 'on');
                    set(obj.goh.port, 'Enable', 'on');
                    set(obj.goh.text_IP_address, 'Enable', 'on');
                    set(obj.goh.text_port, 'Enable', 'on');
                    obj.toggleNTRIP();
                                    
                    %enable approximate position
                    set(obj.goh.text_approx_pos, 'Enable', 'on');
                    set(obj.goh.approx_lat, 'Enable', 'on');
                    set(obj.goh.approx_lon, 'Enable', 'on');
                    set(obj.goh.approx_h, 'Enable', 'on');
                    set(obj.goh.text_approx_lat, 'Enable', 'on');
                    set(obj.goh.text_approx_lon, 'Enable', 'on');
                    set(obj.goh.text_approx_h, 'Enable', 'on');
                    set(obj.goh.text_approx_lat_unit, 'Enable', 'on');
                    set(obj.goh.text_approx_lon_unit, 'Enable', 'on');
                    set(obj.goh.text_approx_h_unit, 'Enable', 'on');
                    
                    set(obj.goh.stopGOstop, 'Enable', 'off');
                    set(obj.goh.text_stopGOstop, 'Enable', 'off');
                    
                elseif (strcmp(contents{get(obj.goh.nav_mon,'Value')},'Rover and Master monitor'))
                    
                    set(obj.goh.text_num_receivers, 'Enable', 'off');
                    set(obj.goh.num_receivers, 'Value', 1);
                    set(obj.goh.num_receivers, 'Enable', 'off');
                    set(obj.goh.com_select_0, 'Enable', 'on');
                    set(obj.goh.com_select_1, 'Enable', 'off');
                    set(obj.goh.com_select_2, 'Enable', 'off');
                    set(obj.goh.com_select_3, 'Enable', 'off');
                    set(obj.goh.protocol_select_0, 'Enable', 'on');
                    set(obj.goh.protocol_select_1, 'Enable', 'off');
                    set(obj.goh.protocol_select_2, 'Enable', 'off');
                    set(obj.goh.protocol_select_3, 'Enable', 'off');
                    set(obj.goh.use_ntrip, 'Enable', 'on');
                    
                    %enable master connection
                    set(obj.goh.IP_address, 'Enable', 'on');
                    set(obj.goh.port, 'Enable', 'on');
                    set(obj.goh.text_IP_address, 'Enable', 'on');
                    set(obj.goh.text_port, 'Enable', 'on');
                    obj.toggleNTRIP();
                                    
                    %disable approximate position
                    set(obj.goh.text_approx_pos, 'Enable', 'off');
                    set(obj.goh.approx_lat, 'Enable', 'off');
                    set(obj.goh.approx_lon, 'Enable', 'off');
                    set(obj.goh.approx_h, 'Enable', 'off');
                    set(obj.goh.text_approx_lat, 'Enable', 'off');
                    set(obj.goh.text_approx_lon, 'Enable', 'off');
                    set(obj.goh.text_approx_h, 'Enable', 'off');
                    set(obj.goh.text_approx_lat_unit, 'Enable', 'off');
                    set(obj.goh.text_approx_lon_unit, 'Enable', 'off');
                    set(obj.goh.text_approx_h_unit, 'Enable', 'off');
                    
                    set(obj.goh.stopGOstop, 'Enable', 'on');
                    set(obj.goh.text_stopGOstop, 'Enable', 'on');
                    obj.toggleStopGoStop();
                end
            end
        end
        
        % This function manage the list box to select the processing type
        function selectProcessingType(obj)            
            contents = cellstr(get(obj.goh.kalman_ls,'String'));
            if (strcmp(contents{get(obj.goh.kalman_ls,'Value')},'Kalman filter'))
                cell_contents = cell(4,1);
                cell_contents{1} = 'Code stand-alone';
                cell_contents{2} = 'Code double difference';
                cell_contents{3} = 'Code and phase stand-alone';
                cell_contents{4} = 'Code and phase double difference';
                cell_contents{5} = 'Code and phase double difference for several receivers';
                set(obj.goh.code_dd_sa, 'String', cell_contents);
                
                set(obj.goh.std_code, 'Enable', 'on');
                set(obj.goh.text_std_code, 'Enable', 'on');
                set(obj.goh.text_std_code_unit, 'Enable', 'on');
                set(obj.goh.std_init, 'Enable', 'on');
                set(obj.goh.text_std_init, 'Enable', 'on');
                set(obj.goh.text_std_init_unit, 'Enable', 'on');
                
                set(obj.goh.toggle_std_phase, 'Enable', 'on');
                set(obj.goh.toggle_std_dtm, 'Enable', 'on');
                obj.togglePhaseStd();
                obj.toggleDtmStd();

                obj.toggleLinearConstraint();
                obj.selectProcessingSolution();
                
                set(obj.goh.cut_off, 'Enable', 'on');
                set(obj.goh.text_cut_off, 'Enable', 'on');
                set(obj.goh.text_cut_off_unit, 'Enable', 'on');
                set(obj.goh.snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres_unit, 'Enable', 'on');
                set(obj.goh.min_sat, 'Enable', 'on');
                set(obj.goh.text_min_sat, 'Enable', 'on');
                set(obj.goh.dyn_mod, 'Enable', 'on');
                %set(obj.goh.text_dyn_mod, 'Enable', 'on');
                
                obj.selectDynMode();
                obj.toggleStopGoStop();
            else
                cell_contents = cell(2,1);
                cell_contents{1} = 'Code stand-alone';
                cell_contents{2} = 'Code double difference';
                old_value = get(obj.goh.code_dd_sa, 'Value');
                if (old_value == 3), set(obj.goh.code_dd_sa, 'Value', 1); end
                if (old_value == 4), set(obj.goh.code_dd_sa, 'Value', 2); end
                set(obj.goh.code_dd_sa, 'String', cell_contents);
                
                obj.selectProcessingSolution();                
                %disable Kalman filters settings
                set(obj.goh.std_X, 'Enable', 'off');
                set(obj.goh.std_Y, 'Enable', 'off');
                set(obj.goh.std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X, 'Enable', 'off');
                set(obj.goh.text_std_Y, 'Enable', 'off');
                set(obj.goh.text_std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X_unit, 'Enable', 'off');
                set(obj.goh.text_std_Y_unit, 'Enable', 'off');
                set(obj.goh.text_std_Z_unit, 'Enable', 'off');
                set(obj.goh.std_code, 'Enable', 'off');
                set(obj.goh.std_phase, 'Enable', 'off');
                set(obj.goh.text_std_code, 'Enable', 'off');
                set(obj.goh.toggle_std_phase, 'Enable', 'off');
                set(obj.goh.text_std_code_unit, 'Enable', 'off');
                set(obj.goh.text_std_phase_unit, 'Enable', 'off');
                set(obj.goh.std_init, 'Enable', 'off');
                set(obj.goh.std_dtm, 'Enable', 'off');
                set(obj.goh.std_vel, 'Enable', 'off');
                set(obj.goh.text_std_init, 'Enable', 'off');
                set(obj.goh.toggle_std_dtm, 'Enable', 'off');
                set(obj.goh.text_std_vel, 'Enable', 'off');
                set(obj.goh.text_std_init_unit, 'Enable', 'off');
                set(obj.goh.text_std_dtm_unit, 'Enable', 'off');
                set(obj.goh.text_std_vel_unit, 'Enable', 'off');
                set(obj.goh.snr_thres, 'Enable', 'off');
                set(obj.goh.text_snr_thres, 'Enable', 'off');
                set(obj.goh.text_snr_thres_unit, 'Enable', 'off');
                set(obj.goh.cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh_unit, 'Enable', 'off');
                set(obj.goh.min_sat, 'Enable', 'off');
                set(obj.goh.text_min_sat, 'Enable', 'off');
                set(obj.goh.antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h_unit, 'Enable', 'off');
                set(obj.goh.dtm_path, 'Enable', 'off');
                set(obj.goh.text_dtm_path, 'Enable', 'off');
                set(obj.goh.browse_dtm_path, 'Enable', 'off');
                set(obj.goh.stopGOstop, 'Enable', 'off');
                set(obj.goh.text_stopGOstop, 'Enable', 'off');
                set(obj.goh.dyn_mod, 'Enable', 'off');
                %set(obj.goh.text_dyn_mod, 'Enable', 'off');
                set(obj.goh.flag_doppler, 'Enable', 'off');
                set(obj.goh.amb_select, 'Enable', 'off');
            end
        end
        
        % This function manage the list box to select the processing kind
        % of solution: Double differences, code + phase, ecc...
        function selectProcessingSolution(obj)            
            contents = cellstr(get(obj.goh.code_dd_sa,'String'));
            % set default values
            set(obj.goh.text_RINEX_rover_obs, 'String', 'RINEX rover observation');
            set(obj.goh.gogps_data, 'Enable', 'on');
            
            if strcmp(contents{get(obj.goh.nav_mon,'Value')},'Code and phase double difference') || strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference for several receivers')
                check_mode = cellstr(get(obj.goh.mode,'String'));
                if (~strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time')) && ...
                        (get(obj.goh.plotproc,'Value'))
                    set(obj.goh.plot_amb, 'Enable', 'on');
                    obj.toggleAmbiguities();
                end
                set(obj.goh.cs_thresh, 'Enable', 'on');
                set(obj.goh.text_cs_thresh, 'Enable', 'on');
                set(obj.goh.text_cs_thresh_unit, 'Enable', 'on');
                set(obj.goh.toggle_std_phase, 'Enable', 'on');
                obj.togglePhaseStd();
                set(obj.goh.snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres_unit, 'Enable', 'on');
                set(obj.goh.constraint, 'Enable', 'on');
                set(obj.goh.ref_path, 'Enable', 'on');
                obj.toggleReferencePath();
                set(obj.goh.stopGOstop, 'Enable', 'on');
                set(obj.goh.text_stopGOstop, 'Enable', 'on');
                obj.toggleStopGoStop();
                cell_contents = cell(4,1);
                cell_contents{1} = 'Const. velocity';
                cell_contents{2} = 'Const. acceleration';
                cell_contents{3} = 'Static';
                cell_contents{4} = 'Variable';
                set(obj.goh.dyn_mod, 'String', cell_contents);
                if (~strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time'))
                    set(obj.goh.flag_doppler, 'Enable', 'on');
                end
                set(obj.goh.amb_select, 'Enable', 'on');
                if strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference for several receivers')
                    %in case of several receivers set RINEX mode (only)
                    %         mode_data = 0;
                    set(obj.goh.rinex_files, 'Value', 1);
                    if(get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files & ~strcmp(contents{get(obj.goh.mode,'Value')},'Real-time'));
                        % enable rover ini file
                        set(obj.goh.RINEX_rover_obs, 'Enable', 'on');
                        set(obj.goh.text_RINEX_rover_obs, 'Enable', 'on');
                        set(obj.goh.text_RINEX_rover_obs, 'String', 'RINEX rover .ini file');
                        set(obj.goh.browse_rover_obs, 'Enable', 'on');
                        % enable master obs file
                        set(obj.goh.RINEX_master_obs, 'Enable', 'on');
                        set(obj.goh.text_RINEX_master_obs, 'Enable', 'on');
                        set(obj.goh.browse_master_obs, 'Enable', 'on');
                        % enable master nav file
                        set(obj.goh.RINEX_nav, 'Enable', 'on');
                        set(obj.goh.text_RINEX_nav, 'Enable', 'on');
                        set(obj.goh.browse_nav, 'Enable', 'on');
                        % disable goGPS binary file
                        set(obj.goh.gogps_data_input, 'Enable', 'off');
                        set(obj.goh.text_gogps_input, 'Enable', 'off');
                        set(obj.goh.browse_gogps_input, 'Enable', 'off');
                    end
                    set(obj.goh.gogps_data, 'Value', 0);
                    set(obj.goh.gogps_data, 'Enable', 'off');
                end
            elseif strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase stand-alone')
                check_mode = cellstr(get(obj.goh.mode,'String'));
                if (~strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time')) & ...
                        (get(obj.goh.plotproc,'Value'))
                    set(obj.goh.plot_amb, 'Enable', 'on');
                    obj.toggleAmbiguities();
                end
                set(obj.goh.cs_thresh, 'Enable', 'on');
                set(obj.goh.text_cs_thresh, 'Enable', 'on');
                set(obj.goh.text_cs_thresh_unit, 'Enable', 'on');
                set(obj.goh.toggle_std_phase, 'Enable', 'on');
                obj.togglePhaseStd();
                set(obj.goh.snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres, 'Enable', 'on');
                set(obj.goh.text_snr_thres_unit, 'Enable', 'on');
                obj.toggleLinearConstraint(0);
                set(obj.goh.constraint, 'Enable', 'off');
                set(obj.goh.stopGOstop, 'Enable', 'off');
                set(obj.goh.text_stopGOstop, 'Enable', 'off');
                set(obj.goh.dyn_mod, 'Enable', 'on');
                cell_contents = cell(3,1);
                cell_contents{1} = 'Const. velocity';
                cell_contents{2} = 'Const. acceleration';
                cell_contents{3} = 'Static';
                old_value = get(obj.goh.dyn_mod, 'Value');
                if (old_value == 4), set(obj.goh.dyn_mod, 'Value', 1); end
                set(obj.goh.dyn_mod, 'String', cell_contents);
                set(obj.goh.flag_doppler, 'Enable', 'on');
                set(obj.goh.amb_select, 'Enable', 'off');
            else
                set(obj.goh.plot_amb, 'Enable', 'off');
                set(obj.goh.no_skyplot_snr, 'Enable', 'on');
                set(obj.goh.plotproc, 'Enable', 'on');
                obj.togglePlotOnProcessing();
                set(obj.goh.cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh, 'Enable', 'off');
                set(obj.goh.text_cs_thresh_unit, 'Enable', 'off');
                set(obj.goh.toggle_std_phase, 'Enable', 'off');
                set(obj.goh.std_phase, 'Enable', 'off');
                set(obj.goh.text_std_phase_unit, 'Enable', 'off');
                obj.toggleLinearConstraint(0);
                set(obj.goh.constraint, 'Enable', 'off');
                set(obj.goh.stopGOstop, 'Enable', 'off');
                set(obj.goh.text_stopGOstop, 'Enable', 'off');
                check_KF = cellstr(get(obj.goh.kalman_ls,'String'));
                if (strcmp(check_KF{get(obj.goh.kalman_ls,'Value')},'Kalman filter'))
                    set(obj.goh.dyn_mod, 'Enable', 'on');
                end
                cell_contents = cell(3,1);
                cell_contents{1} = 'Const. velocity';
                cell_contents{2} = 'Const. acceleration';
                cell_contents{3} = 'Static';
                old_value = get(obj.goh.dyn_mod, 'Value');
                if (old_value == 4), set(obj.goh.dyn_mod, 'Value', 1); end
                set(obj.goh.dyn_mod, 'String', cell_contents);
                set(obj.goh.flag_doppler, 'Enable', 'off');
                set(obj.goh.amb_select, 'Enable', 'off');
            end
            
            if strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase stand-alone') | ...
                    strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code stand-alone')
                
                set(obj.goh.RINEX_master_obs, 'Enable', 'off');
                set(obj.goh.text_RINEX_master_obs, 'Enable', 'off');
                set(obj.goh.browse_master_obs, 'Enable', 'off');
                set(obj.goh.plot_master, 'Enable', 'off');
                set(obj.goh.master_pos, 'Enable', 'off');
                set(obj.goh.crs, 'Enable', 'off');
                set(obj.goh.master_X, 'Enable', 'off');
                set(obj.goh.master_Y, 'Enable', 'off');
                set(obj.goh.master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X, 'Enable', 'off');
                set(obj.goh.text_master_Y, 'Enable', 'off');
                set(obj.goh.text_master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X_unit, 'Enable', 'off');
                set(obj.goh.text_master_Y_unit, 'Enable', 'off');
                set(obj.goh.text_master_Z_unit, 'Enable', 'off');
                set(obj.goh.master_lat, 'Enable', 'off');
                set(obj.goh.master_lon, 'Enable', 'off');
                set(obj.goh.master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat, 'Enable', 'off');
                set(obj.goh.text_master_lon, 'Enable', 'off');
                set(obj.goh.text_master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat_unit, 'Enable', 'off');
                set(obj.goh.text_master_lon_unit, 'Enable', 'off');
                set(obj.goh.text_master_h_unit, 'Enable', 'off');
            else
                contents = cellstr(get(obj.goh.mode,'String'));
                if(get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files & ~strcmp(contents{get(obj.goh.mode,'Value')},'Real-time'));
                    set(obj.goh.RINEX_master_obs, 'Enable', 'on');
                    set(obj.goh.text_RINEX_master_obs, 'Enable', 'on');
                    set(obj.goh.browse_master_obs, 'Enable', 'on');
                end
                set(obj.goh.plot_master, 'Enable', 'on');
                set(obj.goh.master_pos, 'Enable', 'on');
                obj.toggleMasterPosition();
            end
        end
        
        %   INPUT FILE TYPES
        % =================================================================
        
        % Select file type to manage
        function selectFileType(obj, hObject)
            if (strcmp(get(obj.goh.rinex_files, 'Enable'),'on'))
                if (hObject == obj.goh.rinex_files)
                    
                    set(obj.goh.RINEX_rover_obs, 'Enable', 'on');
                    set(obj.goh.RINEX_nav, 'Enable', 'on');
                    set(obj.goh.browse_rover_obs, 'Enable', 'on');
                    set(obj.goh.browse_nav, 'Enable', 'on');
                    set(obj.goh.text_RINEX_rover_obs, 'Enable', 'on');
                    set(obj.goh.text_RINEX_nav, 'Enable', 'on');
                    
                    obj.selectProcessingSolution();
                    
                    set(obj.goh.gogps_data_input, 'Enable', 'off');
                    set(obj.goh.browse_gogps_input, 'Enable', 'off');
                    set(obj.goh.text_gogps_input, 'Enable', 'off');
                    
                    set(obj.goh.stopGOstop, 'Enable', 'off');
                    set(obj.goh.stopGOstop, 'Value', 0);
                    set(obj.goh.text_stopGOstop, 'Enable', 'off');
                    obj.toggleStopGoStop();
                    
                    cell_contents = cell(3,1);
                    cell_contents{1} = 'Const. velocity';
                    cell_contents{2} = 'Const. acceleration';
                    cell_contents{3} = 'Static';
                    old_value = get(obj.goh.dyn_mod, 'Value');
                    if (old_value == 4), set(obj.goh.dyn_mod, 'Value', 1); end
                    set(obj.goh.dyn_mod, 'String', cell_contents);
                    
                else
                    set(obj.goh.RINEX_rover_obs, 'Enable', 'off');
                    set(obj.goh.RINEX_master_obs, 'Enable', 'off');
                    set(obj.goh.RINEX_nav, 'Enable', 'off');
                    set(obj.goh.browse_rover_obs, 'Enable', 'off');
                    set(obj.goh.browse_master_obs, 'Enable', 'off');
                    set(obj.goh.browse_nav, 'Enable', 'off');
                    set(obj.goh.text_RINEX_rover_obs, 'Enable', 'off');
                    set(obj.goh.text_RINEX_master_obs, 'Enable', 'off');
                    set(obj.goh.text_RINEX_nav, 'Enable', 'off');
                    
                    set(obj.goh.gogps_data_input, 'Enable', 'on');
                    set(obj.goh.browse_gogps_input, 'Enable', 'on');
                    set(obj.goh.text_gogps_input, 'Enable', 'on');
                    
                    set(obj.goh.stopGOstop, 'Enable', 'on');
                    set(obj.goh.text_stopGOstop, 'Enable', 'on');
                    obj.toggleStopGoStop();
                    
                    obj.selectProcessingSolution();
                end
            end
        end
        
        %   OPTIONS
        % =================================================================
        
        % Option linear constraint management
        function toggleLinearConstraint(obj, value)
            if nargin == 2
                set(obj.goh.constraint, 'Value', value);
            else
                value = get(obj.goh.constraint,'Value');
            end
            
            if (value)
                set(obj.goh.err_ellipse, 'Enable', 'off');
                set(obj.goh.std_vel, 'Enable', 'on');
                set(obj.goh.text_std_vel, 'Enable', 'on');
                set(obj.goh.text_std_vel_unit, 'Enable', 'on');
            else
                if (get(obj.goh.plotproc,'Value'))
                    set(obj.goh.err_ellipse, 'Enable', 'on');
                end
                set(obj.goh.std_vel, 'Enable', 'off');
                set(obj.goh.text_std_vel, 'Enable', 'off');
                set(obj.goh.text_std_vel_unit, 'Enable', 'off');
            end
        end
        
        % Option Reference Path
        function toggleReferencePath(obj, value)
            if nargin == 2
                set(obj.goh.ref_path, 'Value', value);
            else
                value = get(obj.goh.ref_path,'Value');
            end
            
            contents = cellstr(get(obj.goh.code_dd_sa,'String'));
            if (value) && (strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference') || (strcmp(contents{get(obj.goh.code_dd_sa,'Value')},'Code and phase double difference for several receivers')) )
                set(obj.goh.ref_path_input, 'Enable', 'on');
                set(obj.goh.text_ref_path_input, 'Enable', 'on');
                set(obj.goh.browse_ref_path_input, 'Enable', 'on');
                set(obj.goh.constraint, 'Enable', 'on');
            else
                set(obj.goh.ref_path_input, 'Enable', 'off');
                set(obj.goh.text_ref_path_input, 'Enable', 'off');
                set(obj.goh.browse_ref_path_input, 'Enable', 'off');
                obj.toggleLinearConstraint(0);
                set(obj.goh.constraint, 'Enable', 'off');
            end
        end
        
        % Option Plot while processing
        function togglePlotOnProcessing(obj, value)
            if nargin == 2
                set(obj.goh.plotproc, 'Value', value);
            else
                value = get(obj.goh.plotproc, 'Value');
            end
            
            if (value)
                set(obj.goh.no_skyplot_snr, 'Enable', 'on');
                set(obj.goh.google_earth, 'Enable', 'on');
                set(obj.goh.err_ellipse, 'Enable', 'on');
                set(obj.goh.plot_master, 'Enable', 'on');
                check_mode = cellstr(get(obj.goh.mode,'String'));
                check_phase = cellstr(get(obj.goh.code_dd_sa,'String'));
                if (~strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time') && (strcmp(check_phase{get(obj.goh.code_dd_sa,'Value')}, ...
                        'Code and phase double difference') || strcmp(check_phase{get(obj.goh.code_dd_sa,'Value')},'Code and phase stand-alone')))
                    set(obj.goh.plot_amb, 'Enable', 'on');
                    obj.toggleAmbiguities();
                end
            else
                set(obj.goh.no_skyplot_snr, 'Enable', 'off');
                set(obj.goh.google_earth, 'Enable', 'off');
                set(obj.goh.err_ellipse, 'Enable', 'off');
                set(obj.goh.plot_master, 'Enable', 'off');
                set(obj.goh.plot_amb, 'Enable', 'off');
            end                
        end
        
        % Option Plot ambiguities
        function toggleAmbiguities(obj, value)
            if nargin == 2
                set(obj.goh.plot_amb, 'Value', value);
            else
                value = get(obj.goh.plot_amb,'Value');
            end
            
            if (value)
                obj.enableSkyPlotSNR();
            else
                obj.disableSkyPlotSNR();                
            end
        end      
        
        % Option Use NTRIP
        function toggleNTRIP(obj, value)
            if nargin == 2
                set(obj.goh.use_ntrip, 'Value', value);
            else
                value = get(obj.goh.use_ntrip,'Value');
            end
            
            if (value)
                %enable NTRIP parameters
                set(obj.goh.mountpoint, 'Enable', 'on');
                set(obj.goh.username, 'Enable', 'on');
                set(obj.goh.password, 'Enable', 'on');
                set(obj.goh.show_password, 'Enable', 'on');
                set(obj.goh.text_mountpoint, 'Enable', 'on');
                set(obj.goh.text_username, 'Enable', 'on');
                set(obj.goh.text_password, 'Enable', 'on');
            else
                %disable NTRIP parameters
                set(obj.goh.mountpoint, 'Enable', 'off');
                set(obj.goh.username, 'Enable', 'off');
                set(obj.goh.password, 'Enable', 'off');
                set(obj.goh.show_password, 'Enable', 'off');
                set(obj.goh.text_mountpoint, 'Enable', 'off');
                set(obj.goh.text_username, 'Enable', 'off');
                set(obj.goh.text_password, 'Enable', 'off');
                
                %disable approximate position
                set(obj.goh.text_approx_pos, 'Enable', 'off');
                set(obj.goh.approx_lat, 'Enable', 'off');
                set(obj.goh.approx_lon, 'Enable', 'off');
                set(obj.goh.approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat, 'Enable', 'off');
                set(obj.goh.text_approx_lon, 'Enable', 'off');
                set(obj.goh.text_approx_h, 'Enable', 'off');
                set(obj.goh.text_approx_lat_unit, 'Enable', 'off');
                set(obj.goh.text_approx_lon_unit, 'Enable', 'off');
                set(obj.goh.text_approx_h_unit, 'Enable', 'off');
            end
        end      

        %   KALMANN FILTER
        % =================================================================
        
        % Toggle sto go stop direction estimtion
        function toggleStopGoStop(obj)
            check_nav = cellstr(get(obj.goh.nav_mon,'String'));
            check_KF = cellstr(get(obj.goh.kalman_ls,'String'));
            if (get(obj.goh.stopGOstop,'Value') || strcmp(check_nav{get(obj.goh.nav_mon,'Value')},'Master monitor') || ~strcmp(check_KF{get(obj.goh.kalman_ls,'Value')},'Kalman filter'))
                set(obj.goh.dyn_mod, 'Enable', 'off');
            else
                set(obj.goh.dyn_mod, 'Enable', 'on');
            end
        end
        
        % Select the dynamic mode of the kalman filter
        function selectDynMode(obj)
            global order
            
            contents = cellstr(get(obj.goh.dyn_mod,'String'));
            if (strcmp(contents{get(obj.goh.dyn_mod,'Value')},'Static'))
                order = 1;
                set(obj.goh.std_X, 'Enable', 'off');
                set(obj.goh.std_Y, 'Enable', 'off');
                set(obj.goh.std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X, 'Enable', 'off');
                set(obj.goh.text_std_Y, 'Enable', 'off');
                set(obj.goh.text_std_Z, 'Enable', 'off');
                set(obj.goh.text_std_X_unit, 'Enable', 'off');
                set(obj.goh.text_std_Y_unit, 'Enable', 'off');
                set(obj.goh.text_std_Z_unit, 'Enable', 'off');
            else
                mode = cellstr(get(obj.goh.nav_mon,'String'));
                if (strcmp(mode{get(obj.goh.nav_mon,'Value')},'Navigation'))
                    set(obj.goh.std_X, 'Enable', 'on');
                    set(obj.goh.std_Y, 'Enable', 'on');
                    set(obj.goh.std_Z, 'Enable', 'on');
                    set(obj.goh.text_std_X, 'Enable', 'on');
                    set(obj.goh.text_std_Y, 'Enable', 'on');
                    set(obj.goh.text_std_Z, 'Enable', 'on');
                    set(obj.goh.text_std_X_unit, 'Enable', 'on');
                    set(obj.goh.text_std_Y_unit, 'Enable', 'on');
                    set(obj.goh.text_std_Z_unit, 'Enable', 'on');
                end
                if (strcmp(contents{get(obj.goh.dyn_mod,'Value')},'Const. acceleration'))
                    order = 3;
                elseif (strcmp(contents{get(obj.goh.dyn_mod,'Value')},'Const. velocity'))
                    order = 2;
                else
                    order = 1;
                end
            end
        end
        
        %   KALMANN FILTER - ERROR STD
        % =================================================================

        % Toggle on/off Phase Std parameter       
        function togglePhaseStd(obj, value)
            if nargin == 2
                set(obj.goh.toggle_std_phase, 'Value', value);
            else
                value = get(obj.goh.toggle_std_phase,'Value');
            end

            if (value)
                set(obj.goh.std_phase, 'Enable', 'on');
                set(obj.goh.text_std_phase_unit, 'Enable', 'on');
            else
                set(obj.goh.std_phase, 'Enable', 'off');
                set(obj.goh.text_std_phase_unit, 'Enable', 'off');
            end
        end
        
        % Toggle on/off DTM Std parameter
        function toggleDtmStd(obj, value)
            if nargin == 2
                set(obj.goh.toggle_std_dtm, 'Value', value);
            else
                value = get(obj.goh.toggle_std_dtm,'Value');
            end

            if (value)
                set(obj.goh.std_dtm, 'Enable', 'on');
                set(obj.goh.text_std_dtm_unit, 'Enable', 'on');
                set(obj.goh.antenna_h, 'Enable', 'on');
                set(obj.goh.text_antenna_h, 'Enable', 'on');
                set(obj.goh.text_antenna_h_unit, 'Enable', 'on');
                set(obj.goh.dtm_path, 'Enable', 'on');
                set(obj.goh.text_dtm_path, 'Enable', 'on');
                set(obj.goh.browse_dtm_path, 'Enable', 'on');
            else
                set(obj.goh.std_dtm, 'Enable', 'off');
                set(obj.goh.text_std_dtm_unit, 'Enable', 'off');
                set(obj.goh.antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h, 'Enable', 'off');
                set(obj.goh.text_antenna_h_unit, 'Enable', 'off');
                set(obj.goh.dtm_path, 'Enable', 'off');
                set(obj.goh.text_dtm_path, 'Enable', 'off');
                set(obj.goh.browse_dtm_path, 'Enable', 'off');
            end
        end
        
        %   KALMANN FILTER - AMBIGUITY RESTART METHOD
        % =================================================================
        
        % Select the method to be used to fix the ambiguity after anomalies
        function selectAmbiguityRestartMethod(obj)
            global amb_restart_method
            contents = cellstr(get(obj.goh.amb_select,'String'));
            selection = contents{get(obj.goh.amb_select,'Value')};
            if (strcmp(selection, 'Observed code - phase difference'))
                amb_restart_method = 0;
            elseif (strcmp(selection, 'Kalman-predicted code - phase difference'))
                amb_restart_method = 1;
            else
                amb_restart_method = 2;
            end
        end        
        
        %   MASTER STATION
        % =================================================================
            
        % Master station source: RINEX position / interface
        function toggleMasterPosition(obj, value)
            if nargin == 2
                set(obj.goh.master_pos, 'Value', value);
            else
                value = get(obj.goh.master_pos,'Value');
            end
            
            if (value)
                set(obj.goh.crs, 'Enable', 'off');
                set(obj.goh.master_X, 'Enable', 'off');
                set(obj.goh.master_Y, 'Enable', 'off');
                set(obj.goh.master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X, 'Enable', 'off');
                set(obj.goh.text_master_Y, 'Enable', 'off');
                set(obj.goh.text_master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X_unit, 'Enable', 'off');
                set(obj.goh.text_master_Y_unit, 'Enable', 'off');
                set(obj.goh.text_master_Z_unit, 'Enable', 'off');
                set(obj.goh.master_lat, 'Enable', 'off');
                set(obj.goh.master_lon, 'Enable', 'off');
                set(obj.goh.master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat, 'Enable', 'off');
                set(obj.goh.text_master_lon, 'Enable', 'off');
                set(obj.goh.text_master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat_unit, 'Enable', 'off');
                set(obj.goh.text_master_lon_unit, 'Enable', 'off');
                set(obj.goh.text_master_h_unit, 'Enable', 'off');
            else
                set(obj.goh.crs, 'Enable', 'on');
                obj.toggleMasterPositionType();
            end
        end
        
        % Enable / disable the right coordinate type to read
        function toggleMasterPositionType(obj)
            contents = cellstr(get(obj.goh.crs,'String'));
            if (strcmp(contents{get(obj.goh.crs,'Value')},'ECEF (X,Y,Z)'))
                set(obj.goh.master_X, 'Enable', 'on');
                set(obj.goh.master_Y, 'Enable', 'on');
                set(obj.goh.master_Z, 'Enable', 'on');
                set(obj.goh.text_master_X, 'Enable', 'on');
                set(obj.goh.text_master_Y, 'Enable', 'on');
                set(obj.goh.text_master_Z, 'Enable', 'on');
                set(obj.goh.text_master_X_unit, 'Enable', 'on');
                set(obj.goh.text_master_Y_unit, 'Enable', 'on');
                set(obj.goh.text_master_Z_unit, 'Enable', 'on');
                set(obj.goh.master_lat, 'Enable', 'off');
                set(obj.goh.master_lon, 'Enable', 'off');
                set(obj.goh.master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat, 'Enable', 'off');
                set(obj.goh.text_master_lon, 'Enable', 'off');
                set(obj.goh.text_master_h, 'Enable', 'off');
                set(obj.goh.text_master_lat_unit, 'Enable', 'off');
                set(obj.goh.text_master_lon_unit, 'Enable', 'off');
                set(obj.goh.text_master_h_unit, 'Enable', 'off');
            else
                set(obj.goh.master_X, 'Enable', 'off');
                set(obj.goh.master_Y, 'Enable', 'off');
                set(obj.goh.master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X, 'Enable', 'off');
                set(obj.goh.text_master_Y, 'Enable', 'off');
                set(obj.goh.text_master_Z, 'Enable', 'off');
                set(obj.goh.text_master_X_unit, 'Enable', 'off');
                set(obj.goh.text_master_Y_unit, 'Enable', 'off');
                set(obj.goh.text_master_Z_unit, 'Enable', 'off');
                set(obj.goh.master_lat, 'Enable', 'on');
                set(obj.goh.master_lon, 'Enable', 'on');
                set(obj.goh.master_h, 'Enable', 'on');
                set(obj.goh.text_master_lat, 'Enable', 'on');
                set(obj.goh.text_master_lon, 'Enable', 'on');
                set(obj.goh.text_master_h, 'Enable', 'on');
                set(obj.goh.text_master_lat_unit, 'Enable', 'on');
                set(obj.goh.text_master_lon_unit, 'Enable', 'on');
                set(obj.goh.text_master_h_unit, 'Enable', 'on');
            end
        end

        %   ROVER PORT AND PROTOCOL
        % =================================================================
        
        % Set number of phisical receiver to connect with
        function setNumReceiverIn(obj)
            contents = cellstr(get(obj.goh.num_receivers,'String'));
            check_mode = cellstr(get(obj.goh.mode,'String'));
            contents_nav_mon = cellstr(get(obj.goh.nav_mon,'String'));
            if (strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time') && strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Rover monitor'))
                set(obj.goh.text_num_receivers, 'Enable', 'on');
                set(obj.goh.num_receivers, 'Enable', 'on');
                if (size(contents,1) >= get(obj.goh.num_receivers,'Value'))
                    if (strcmp(contents{get(obj.goh.num_receivers,'Value')},'1'))
                        set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                        set(obj.goh.com_select_1, 'Enable', 'off'); set(obj.goh.protocol_select_1, 'Enable', 'off');
                        set(obj.goh.com_select_2, 'Enable', 'off'); set(obj.goh.protocol_select_2, 'Enable', 'off');
                        set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
                    elseif (strcmp(contents{get(obj.goh.num_receivers,'Value')},'2'))
                        set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                        set(obj.goh.com_select_1, 'Enable', 'on');  set(obj.goh.protocol_select_1, 'Enable', 'on');
                        set(obj.goh.com_select_2, 'Enable', 'off'); set(obj.goh.protocol_select_2, 'Enable', 'off');
                        set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
                    elseif (strcmp(contents{get(obj.goh.num_receivers,'Value')},'3'))
                        set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                        set(obj.goh.com_select_1, 'Enable', 'on');  set(obj.goh.protocol_select_1, 'Enable', 'on');
                        set(obj.goh.com_select_2, 'Enable', 'on');  set(obj.goh.protocol_select_2, 'Enable', 'on');
                        set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
                    elseif (strcmp(contents{get(obj.goh.num_receivers,'Value')},'4'))
                        set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                        set(obj.goh.com_select_1, 'Enable', 'on');  set(obj.goh.protocol_select_1, 'Enable', 'on');
                        set(obj.goh.com_select_2, 'Enable', 'on');  set(obj.goh.protocol_select_2, 'Enable', 'on');
                        set(obj.goh.com_select_3, 'Enable', 'on');  set(obj.goh.protocol_select_3, 'Enable', 'on');
                    end
                else
                    set(obj.goh.num_receivers,'Value',1);
                    set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                    set(obj.goh.com_select_1, 'Enable', 'off'); set(obj.goh.protocol_select_1, 'Enable', 'off');
                    set(obj.goh.com_select_2, 'Enable', 'off'); set(obj.goh.protocol_select_2, 'Enable', 'off');
                    set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
                end
            elseif (strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time') && strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Master monitor'))
                set(obj.goh.text_num_receivers, 'Enable', 'off');
                set(obj.goh.num_receivers, 'Enable', 'off');
                set(obj.goh.com_select_0, 'Enable', 'off');
                set(obj.goh.com_select_1, 'Enable', 'off');
                set(obj.goh.com_select_2, 'Enable', 'off');
                set(obj.goh.com_select_3, 'Enable', 'off');
                set(obj.goh.protocol_select_0, 'Enable', 'off');
                set(obj.goh.protocol_select_1, 'Enable', 'off');
                set(obj.goh.protocol_select_2, 'Enable', 'off');
                set(obj.goh.protocol_select_3, 'Enable', 'off');
            elseif (strcmp(check_mode{get(obj.goh.mode,'Value')},'Real-time') && (strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Navigation') | strcmp(contents_nav_mon{get(obj.goh.nav_mon,'Value')},'Rover and Master monitor')))
                set(obj.goh.num_receivers,'Value',1);
                set(obj.goh.com_select_0, 'Enable', 'on');  set(obj.goh.protocol_select_0, 'Enable', 'on');
                set(obj.goh.com_select_1, 'Enable', 'off'); set(obj.goh.protocol_select_1, 'Enable', 'off');
                set(obj.goh.com_select_2, 'Enable', 'off'); set(obj.goh.protocol_select_2, 'Enable', 'off');
                set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
            else
                set(obj.goh.num_receivers,'Value',1);
                set(obj.goh.com_select_0, 'Enable', 'off');  set(obj.goh.protocol_select_0, 'Enable', 'off');
                set(obj.goh.com_select_1, 'Enable', 'off'); set(obj.goh.protocol_select_1, 'Enable', 'off');
                set(obj.goh.com_select_2, 'Enable', 'off'); set(obj.goh.protocol_select_2, 'Enable', 'off');
                set(obj.goh.com_select_3, 'Enable', 'off'); set(obj.goh.protocol_select_3, 'Enable', 'off');
            end
        end
        
        %   SKY PLOT SNR
        % =================================================================
        
        function enableSkyPlotSNR(obj)
            set(obj.goh.no_skyplot_snr, 'Enable', 'on');
        end        
        function disableSkyPlotSNR(obj)
            set(obj.goh.no_skyplot_snr, 'Enable', 'off');
        end

        
        
        
        
                % Load the state of the gui from file.
        function loadState(obj,filename)            
            load(filename); % the file contains the variable state
            obj.status = state;
            
            set(obj.goh.master_pos, 'Value', state.master_pos);
            set(obj.goh.constraint, 'Value', state.constraint);
            set(obj.goh.ref_path, 'Value', state.ref_path);
            set(obj.goh.plot_master, 'Value', state.plot_master);
            set(obj.goh.google_earth, 'Value', state.google_earth);
            set(obj.goh.err_ellipse, 'Value', state.err_ellipse);
            set(obj.goh.use_ntrip, 'Value', state.use_ntrip);
            set(obj.goh.plot_amb, 'Value', state.plot_amb);
            set(obj.goh.no_skyplot_snr, 'Value', state.no_skyplot_snr);
            set(obj.goh.plotproc, 'Value', state.plotproc);
            set(obj.goh.RINEX_rover_obs,'String', state.RINEX_rover_obs);
            set(obj.goh.RINEX_master_obs,'String', state.RINEX_master_obs);
            set(obj.goh.RINEX_nav,'String', state.RINEX_nav);
            set(obj.goh.gogps_data_input,'String', state.gogps_data_input);
            set(obj.goh.gogps_data_output,'String', state.gogps_data_output);
            set(obj.goh.gogps_data_output_prefix,'String', state.gogps_data_output_prefix);
            set(obj.goh.ref_path_input,'String', state.ref_path_input);
            set(obj.goh.dtm_path,'String', state.dtm_path);
            set(obj.goh.master_X,'String', state.master_X);
            set(obj.goh.master_Y,'String', state.master_Y);
            set(obj.goh.master_Z,'String', state.master_Z);
            set(obj.goh.master_lat,'String', state.master_lat);
            set(obj.goh.master_lon,'String', state.master_lon);
            set(obj.goh.master_h,'String', state.master_h);
            set(obj.goh.crs,'Value', state.crs);
            set(obj.goh.mode,'Value', state.mode);
            set(obj.goh.nav_mon,'Value', state.nav_mon);
            set(obj.goh.kalman_ls,'Value', state.kalman_ls);
            set(obj.goh.code_dd_sa,'Value', state.code_dd_sa);
            set(obj.goh.rinex_files,'Value', state.rinex_files);
            set(obj.goh.gogps_data,'Value', state.gogps_data);
            set(obj.goh.std_X,'String', state.std_X);
            set(obj.goh.std_Y,'String', state.std_Y);
            set(obj.goh.std_Z,'String', state.std_Z);
            set(obj.goh.std_code,'String', state.std_code);
            set(obj.goh.std_phase,'String', state.std_phase);
            set(obj.goh.std_dtm,'String', state.std_dtm);
            set(obj.goh.toggle_std_phase,'Value', state.toggle_std_phase);
            set(obj.goh.toggle_std_dtm,'Value', state.toggle_std_dtm);
            set(obj.goh.std_init,'String', state.std_init);
            set(obj.goh.std_vel,'String', state.std_vel);
            set(obj.goh.cs_thresh,'String', state.cs_thresh);
            set(obj.goh.flag_doppler,'Value', state.flag_doppler);
            set(obj.goh.amb_select,'Value', state.amb_select);
            set(obj.goh.cut_off,'String', state.cut_off);
            set(obj.goh.snr_thres,'String', state.snr_thres);
            set(obj.goh.antenna_h,'String', state.antenna_h);
            set(obj.goh.min_sat,'String', state.min_sat);
            set(obj.goh.dyn_mod,'Value', state.dyn_mod);
            set(obj.goh.weight_0,'Value', state.weight_0);
            set(obj.goh.weight_1,'Value', state.weight_1);
            set(obj.goh.weight_2,'Value', state.weight_2);
            set(obj.goh.weight_3,'Value', state.weight_3);
            contents = get(obj.goh.com_select_0,'String');
            select_0 = 1; select_1 = 1; select_2 = 1; select_3 = 1;
            for i = 1 : numel(contents)
                if (strcmp(contents{i},state.com_select_0))
                    select_0 = i;
                elseif (strcmp(contents{i},state.com_select_1))
                    select_1 = i;
                elseif (strcmp(contents{i},state.com_select_2))
                    select_2 = i;
                elseif (strcmp(contents{i},state.com_select_3))
                    select_3 = i;
                end
            end
            set(obj.goh.com_select_0,'Value', select_0);
            set(obj.goh.com_select_1,'Value', select_1);
            set(obj.goh.com_select_2,'Value', select_2);
            set(obj.goh.com_select_3,'Value', select_3);
            set(obj.goh.protocol_select_0,'Value', state.protocol_select_0);
            set(obj.goh.protocol_select_1,'Value', state.protocol_select_1);
            set(obj.goh.protocol_select_2,'Value', state.protocol_select_2);
            set(obj.goh.protocol_select_3,'Value', state.protocol_select_3);
            set(obj.goh.num_receivers,'Value', state.num_receivers);
            set(obj.goh.IP_address,'String', state.IP_address);
            set(obj.goh.port,'String', state.port);
            set(obj.goh.mountpoint,'String', state.mountpoint);
            set(obj.goh.username,'String', state.username);
            obj.setPassword(state.password);            
            set(obj.goh.approx_lat,'String', state.approx_lat);
            set(obj.goh.approx_lon,'String', state.approx_lon);
            set(obj.goh.approx_h,'String', state.approx_h);
            set(obj.goh.stopGOstop,'Value', state.stopGOstop);
            
            obj.toggleAmbiguities();
            obj.toggleLinearConstraint();
            obj.toggleMasterPosition();
            obj.selectProcessingType();
            obj.togglePlotOnProcessing();
            obj.selectDynMode();
            obj.selectProcassingApproach();
            obj.setNumReceiverIn();
            obj.toggleStopGoStop();
            if(get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files);
                obj.selectFileType(obj.goh.rinex_files);
            else
                obj.selectFileType(obj.goh.gogps_data);
            end
        end
    end
    
    % =========================================================================
    %    PRIVATE METHODS
    % =========================================================================

    methods(Access = 'private')
        
        function init(obj, handles)
            obj.goh = handles;  % Save the handle of the figure
            
            % Choose default command line output for gui_goGPS_unix
            obj.goh.output = obj.goh.main_panel;
            
            set(obj.goh.main_panel,'CloseRequestFcn',@obj.closeGUI);
            
            % Update handles structure
            guidata(obj.goh.main_panel, obj.goh);
            
            if exist('../data/settings/last_settings.mat','file')
                obj.loadState('../data/settings/last_settings.mat');
            else
                obj.loadState('../data/settings/default_settings.mat');
            end
            
            %pixels
            set(obj.goh.main_panel, 'Units', 'pixels' );
            
            %get display size
            screenSize = get(0, 'ScreenSize');
            
            %calculate the center of the display
            position = get(obj.goh.main_panel, 'Position');
            position(1) = (screenSize(3)-position(3))/2;
            position(2) = (screenSize(4)-position(4))/2;
            
            %center the window
            set(obj.goh.main_panel, 'Position', position);            
        end

        function closeGUI(obj,src,evnt)
            selection = questdlg('Do you want to quit goGPS?',...
                'Close Request Function',...
                'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    delete(gcf)
                case 'No'
                    return
            end
        end

    end
end