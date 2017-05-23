%   CLASS File_Wizard
% =========================================================================
%
% DESCRIPTION
%   Class to check and prepare the files needed by the processing
%   (e.g. navigational files)
%
% EXAMPLE
%   settings = FTP_Server();
%
% FOR A LIST OF CONSTANTs and METHODS use doc File_Wizard
%
% REQUIRES:
%   goGPS settings;
%
% COMMENTS
%   Server structure:

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef File_Wizard < handle

    properties (SetAccess = protected, GetAccess = public)
        % custom ftp downloader
        ftpd_custom;
        % download sources
        source = struct('igscb', struct('ftpd', FTP_Downloader('ftp://igscb.jpl.nasa.gov/','21'), ...
                                        'par',  struct('gps', struct('path', '/pub/gps/', ...
                                                                     'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
                                                                                                    'eph', struct('final', '${WWWW}/igs${WWWWD}.sp3', ...
                                                                                                                  'rapid', '${WWWW}/igr${WWWWD}.sp3', ...
                                                                                                                  'ultra', '${WWWW}/igu${WWWWD}_${6H}.sp3'), ...
                                                                                                    'clk', struct('final', '${WWWW}/igs${WWWWD}.clk', ...
                                                                                                                  'rapid', '${WWWW}/igr${WWWWD}.clk'), ...
                                                                                                    'clk_30s', struct('final', '${WWWW}/igs${WWWWD}.clk_30s'), ...
                                                                                                    'erp', struct('final', '${WWWW}/igs${WWWW}7.erp')))), ...
                                                       'glo', struct('path', '/pub/glonass/products/', ...
                                                                     'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
                                                                                                    'eph', struct('final', '${WWWW}/igl${WWWWD}.sp3', ...
                                                                                                                  'rapid', '${WWWW}/igv${WWWWD}.sp3')))))), ...
                        'cddis', struct('ftpd', FTP_Downloader('ftp://cddis.gsfc.nasa.gov/','21'), ...
                                        'par',  struct('gps', struct('path', '/pub/gps/', ...
                                                                     'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/igs${WWWWD}.sp3', ...
                                                                                                                  'rapid', 'products/${WWWW}/igr${WWWWD}.sp3', ...
                                                                                                                  'ultra', 'products/${WWWW}/igu${WWWWD}_${6H}.sp3', ...
                                                                                                                  'broadcast',  'data/daily/${YYYY}/brdc/brdc${DOY}0.${YY}n'), ...
                                                                                                    'clk', struct('final', 'products/${WWWW}/igs${WWWWD}.clk', ...
                                                                                                                  'rapid', 'products/${WWWW}/igr${WWWWD}.clk'), ...
                                                                                                    'clk_30s', struct('final', 'products/${WWWW}/igs${WWWWD}.clk_30s'), ...
                                                                                                    'erp', struct('final', 'products/${WWWW}/igs${WWWW}7.erp')),...
                                                                                      'cod', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/cod${WWWWD}.eph'), ...
                                                                                                    'clk_05s', struct('final', 'products/${WWWW}/cod${WWWWD}.clk_05s')), ...
                                                                                      'jpl', struct('name', 'Jet Propulsion Laboratory (JPL) - NASA (USA)', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/jpl${WWWWD}.sp3'), ...
                                                                                                    'clk_30s', struct('final', 'products/${WWWW}/jpl${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'products/${WWWW}/jpl${WWWW}7.erp')), ...
                                                                                      'esa', struct('name', 'European Space Operations Centre (ESA) - ESA (Europe)', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/esa${WWWWD}.sp3'), ...
                                                                                                    'clk_30s', struct('final', 'products/${WWWW}/esa${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'products/${WWWW}/esa${WWWW}7.erp')), ...
                                                                                      'gbm', struct('name', 'GeoForschungsZentrum Potsdam (GBM) - GFZ (Germany)', ...
                                                                                                    'eph', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
                                                                                                                  'broadcast',  'data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
                                                                                                    'clk_30s', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'products/mgex/${WWWW}/gbm${WWWW}7.erp'), ...
                                                                                                    'bds', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.bias')) ...
                                                                                      )), ...
                                                       'glo', struct('path', '/', ...
                                                                     'center', struct('igs', struct('name', 'International GNSS System - IGS orbit combination', ...
                                                                                                    'eph', struct('final', 'glonass/products/${WWWW}/igl${WWWWD}.sp3', ...
                                                                                                                  'rapid', 'glonass/products/${WWWW}/igv${WWWWD}.sp3', ...
                                                                                                                  'broadcast',  'pub/gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p')), ...
                                                                                      'emx', struct('name', 'Natural Resources Canada - NRCAN GNSS solution', ...
                                                                                                    'eph', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.sp3'), ...
                                                                                                    'clk_30s', struct('final', 'glonass/products/${WWWW}/emx${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'glonass/products/${WWWW}/emx${WWWW}7.erp')), ...
                                                                                      'cod', struct('name', 'Center for Orbit Determination in Europe (CODE) - AIUB (Switzerland)', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/cod${WWWWD}.eph'), ...
                                                                                                    'clk_05s', struct('final', 'products/${WWWW}/cod${WWWWD}.clk_05s')), ...
                                                                                      'esa', struct('name', 'European Space Operations Centre (ESA) - ESA (Europe)', ...
                                                                                                    'eph', struct('final', 'products/${WWWW}/esa${WWWWD}.sp3'), ...
                                                                                                    'clk_30s', struct('final', 'products/${WWWW}/esa${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'products/${WWWW}/esa${WWWW}7.erp')), ...
                                                                                      'gbm', struct('name', 'GeoForschungsZentrum Potsdam (GFZ)', ...
                                                                                                    'eph', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
                                                                                                                  'broadcast',  '/pub/gnss/data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
                                                                                                    'clk_30s', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWW}7.erp'), ...
                                                                                                    'bds', struct('final', '/pub/gnss/products/mgex/${WWWW}/gbm${WWWWD}.bias')))), ...
                                                       'mxd', struct('path', '/pub/gnss/', ...
                                                                     'center', struct('gbm', struct('name', 'GeoForschungsZentrum Potsdam (GFZ)', ...
                                                                                                    'eph', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.sp3', ...
                                                                                                                  'broadcast',  'data/campaign/mgex/daily/rinex3/${YYYY}/brdm/brdm${DOY}0.${YY}p'), ...
                                                                                                    'clk_30s', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.clk'), ...
                                                                                                    'erp', struct('final', 'products/mgex/${WWWW}/gbm${WWWW}7.erp'), ...
                                                                                                    'bds', struct('final', 'products/mgex/${WWWW}/gbm${WWWWD}.bias')))))));
    	date_start; % first epoch of common observations (among all the obs files)
        date_stop;  % last epoch of common observations (among all the obs files)
    end

    properties (SetAccess = protected, GetAccess = protected)
        state = Go_State.getCurrentSettings();     %  Global state, to import custom server and service preferences
    end

    properties (SetAccess = private, GetAccess = private)
        logger = Logger.getInstance(); % Handler to the logger object
        ftp_downloader;
    end

    methods
        function this = File_Wizard(state)
            % Constructor
            %  SYNTAX File_Wizard(<state>)
            % Uses state for getting settings info
            % Modify state to update eph_name and clk_name

            if (nargin >= 1)
                this.state = handle(state);
            else
                this.state = Go_State.getCurrentSettings();
            end
        end

        function conjureFiles(this, date_start, date_stop)
            % Prepare all the files needed for processing
            if (nargin == 1)
                [date_start, date_stop] = this.conjureObsFile();
            end
            this.state.setProcessingTime(date_start, date_stop, false);
            this.state.updateNavFileName();
            this.conjureNavFiles(date_start, date_stop);
        end

        function [first_epoch, last_epoch] = conjureObsFile(this)
            % Prepare the extended file name of the files to be used in goGPS
            % In a future here I'll download the required navigational files of a station in a network

            first_target_files = this.state.getTargetPath(1);
            fh = File_Rinex(first_target_files);
            this.logger.newLine();
            first_epoch = fh.first_epoch.first;
            last_epoch = fh.last_epoch.last;
        end

        function conjureNavFiles(this, date_start, date_stop)
            % prepare the navigational files needed for processing

            eph_ok = this.state.checkNavEphFiles();
            nav_ok = eph_ok && this.state.checkNavClkFiles();
            if (~nav_ok)
                fnp = File_Name_Processor();

                % Just in case I need it, import custom server
                [addr, port, path, nav_name, clk_name, ~] = this.state.getCustomArchive();
                this.ftpd_custom = FTP_Downloader(addr, port, path);
                clear addr port path;

                % Get lists of preferences
                archive_list = this.state.getNavArchive();
                eph_type_list = this.state.getNavEphType();
                clk_type_list = this.state.getNavClkType();

                active_ss = this.state.cc.getActive();
                if sum(active_ss(3:end))
                    % Multiconstellation orbits are needed
                    provider_list = this.state.getNavMixedProvider();
                    ss = 'mxd';
                elseif (this.state.cc.isGloActive())
                    provider_list = this.state.getNavGloProvider();
                    ss = 'glo';
                else
                    provider_list = this.state.getNavGpsProvider();
                    ss = 'gps';
                end
                clear active_ss;

                % Start the search for navigational files
                % -> stop when found / no more places where to search for them

                t = 0; % type of wanted ephemeris
                while (~nav_ok && (t < numel(eph_type_list)))
                    t = t + 1;
                    eph_type = eph_type_list{t};
                    p = 0; % provider of wanted ephemeris
                    while (~nav_ok && (p < numel(provider_list)))
                        p = p + 1;
                        provider = provider_list{p};
                        a = 0; % archive of wanted ephemeris
                        while (~nav_ok && (a < numel(archive_list)))

                            a = a + 1;
                            archive = archive_list{a};
                            if strcmp(archive, 'custom')
                                % custom provider is selected
                                % download ephemeris
                                step_sec = fnp.getStepSec(nav_name);
                                tmp_date_start = date_start.getCopy; tmp_date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                                tmp_date_stop = date_stop.getCopy; tmp_date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin

                                file_list = flipud(fnp.dateKeyRepBatch(nav_name, tmp_date_start, tmp_date_stop));
                                this.ftpd_custom.download(file_list, this.state.getNavEphDir());
                                [~, name, ext] = fileparts(nav_name);
                                this.state.setNavEphFile(strcat(name, ext));

                                % download clocks
                                step_sec = fnp.getStepSec(clk_name);
                                tmp_date_start = date_start.getCopy; tmp_date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                                tmp_date_stop = date_stop.getCopy; tmp_date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
                                file_list = flipud(fnp.dateKeyRepBatch(clk_name, tmp_date_start, tmp_date_stop));
                                this.ftpd_custom.download(file_list, this.state.getNavClkDir());
                                [~, name, ext] = fileparts(clk_name);
                                this.state.setNavClkFile(strcat(name, ext));

                                % match the name of the ephemeris to use with what I've just downloaded
                                this.state.updateNavFileName();

                                eph_ok = this.state.checkNavEphFiles();
                                nav_ok =  eph_ok && this.state.checkNavClkFiles();
                            elseif isfield(this.source, archive) &&  ...
                                   isfield(this.source.(archive).par, ss) && ...
                                   isfield(this.source.(archive).par.(ss).center, provider) && ...
                                   isfield(this.source.(archive).par.(ss).center.(provider).eph, eph_type)

                                % Download navigational
                                nav_name = this.source.(archive).par.(ss).center.(provider).eph.(eph_type);
                                step_sec = fnp.getStepSec(nav_name);
                                tmp_date_start = date_start.getCopy; tmp_date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                                tmp_date_stop = date_stop.getCopy; tmp_date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
                                file_list = flipud(fnp.dateKeyRepBatch(nav_name, tmp_date_start, tmp_date_stop));
                                this.source.(archive).ftpd.download(this.source.(archive).par.(ss).path, file_list, this.state.getNavEphDir());
                                [~, name, ext] = fileparts(nav_name);
                                this.state.setNavEphFile(strcat(name, ext));

                                % match the name of the ephemeris to use with what I've just downloaded
                                this.state.updateNavFileName();

                                eph_ok = this.state.checkNavEphFiles();

                                % if nav_ok try to download clocks
                                if (eph_ok && (~strcmp(eph_type,'ultra') && ~strcmp(eph_type,'broadcast')))
                                    c = 0;
                                    clk_ok = this.state.checkNavClkFiles();
                                    while (c < numel(clk_type_list)) && ~clk_ok
                                        c = c + 1;
                                        clk_type = clk_type_list{c};
                                        if isfield(this.source.(archive).par.(ss).center.(provider), clk_type) && ...
                                           isfield(this.source.(archive).par.(ss).center.(provider).(clk_type), eph_type)
                                            % Download navigational
                                            clk_name = this.source.(archive).par.(ss).center.(provider).(clk_type).(eph_type);
                                            step_sec = fnp.getStepSec(clk_name);
                                            tmp_date_start = date_start.getCopy; tmp_date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                                            tmp_date_stop = date_stop.getCopy; tmp_date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
                                            file_list = flipud(fnp.dateKeyRepBatch(clk_name, tmp_date_start, tmp_date_stop));
                                            this.source.(archive).ftpd.download(this.source.(archive).par.(ss).path, file_list, this.state.getNavClkDir());
                                            [~, name, ext] = fileparts(clk_name);
                                            this.state.setNavClkFile(strcat(name, ext));

                                            % match the name of the ephemeris to use with what I've just downloaded
                                            this.state.updateNavFileName();

                                            clk_ok = this.state.checkNavClkFiles();
                                        end
                                    end
                                else
                                    clk_ok = true;
                                end
                                if ~clk_ok
                                    clk_ok = this.state.checkNavClkFiles();
                                end
                                nav_ok = eph_ok && clk_ok;
                            end
                        end
                    end
                end
            else
                this.logger.addStatusOk('Navigational files are present ^_^');
                this.logger.newLine();
            end
        end
    end

end
