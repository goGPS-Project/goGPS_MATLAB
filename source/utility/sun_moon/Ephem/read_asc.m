function [ scan, error ] = read_asc( denum, fn, dontSave )
%READ_ASC read an ascii JPL ephemeris file from hard-drive or internet
persistent m;
if isempty(m)
    m = containers.Map('KeyType','char','ValueType','any');
end
if nargin < 3
    dontSave = false;
end

log = Core.getLogger();
error = 1;
if exist(fn, 'file') == 2
    % found locally
    error = 0;
    if m.isKey(fn)
        scan = m(fn);
    else
        log.addMessage(log.indent(sprintf('Reading %s...', fn)));
        EPHFILE = fopen (fn, 'r');
        if EPHFILE ~= -1
            scan = textscan(EPHFILE,'%f');
            if isempty(scan)
                error = 1;
            else
                fclose (EPHFILE);
                scan = scan{1};
                if length(scan) < 4
                    log.addWarning('read_asc: Bad textscan of %s\n',fn);
                    error = 9;
                else
                    if dontSave == false
                        m(fn) = scan;
                    end
                end
            end
        else
            scan = zeros(0,0);
            log.addWarning('read_asc: Could not open %s\n',fn);
            error = 10;
        end
    end
end
if error > 0
    if isnumeric(denum)
        denum = sprintf('%d',denum);
    end
    % retrieve from NASA
    url = sprintf('%s/de%s/%s',Ephem.JPLftpSite,denum,fn);
    if m.isKey(url)
        scan = m(url);
    else
        [goGPS_path] = which('goGPS');
        [goGPS_dir] = fileparts(goGPS_path);
        log.addMessage(log.indent(sprintf('Reading %s...', [goGPS_dir '/reserved/JPL/' fn])));
        % If someone needs to run goGPS from offline for the first time
        % it can download ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de436/ascp01950.436 
        % and ascp02050.436, and put them into "./reserved/JPL/"
        fid = fopen([goGPS_dir '/reserved/JPL/' fn], 'r'); 
        if fid > 0            
            EPHFILE = fread(fid, '*char')';
            status = 1;
            fclose(fid);
        else
            log.addMessage(log.indent(sprintf('File not found, downloading it,\nDownloading "%s"...', url)));
            [EPHFILE, status] = urlread(url);
        end
        if status == 1 && ~isempty(EPHFILE)
            scan = textscan(EPHFILE,'%f');
            scan = scan{1};
            if isempty(scan)
                fprintf(2,'read_asc: Bad textscan of %s\n',url);
                error = 9;
            else
                if dontSave == false
                    m(url) = scan;
                end
            end
        else
            error('Whitout %s sun and moon ephemerides cannot be computed', url);
            scan = zeros(0,0);
            fprintf(2,'read_asc: Could not download %s\n',url);
            error = 10;
        end
    end
end
end

