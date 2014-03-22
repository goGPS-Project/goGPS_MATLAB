function write_RINEX_obs(filename, receiver, antenna, pr1_R, pr2_R, ph1_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, time_R, date, pos_R, interval, flag_P1)

%----------------------------------------------------------------------------------------------
% RINEX OBSERVATION FILE
%----------------------------------------------------------------------------------------------

%displaying
fprintf(['Writing: ' filename '... ']);

%create RINEX observation file
fid_obs = fopen(filename,'wt');

%maximum filename length for MARKER NAME field
pos = find(filename == '/');
if (isempty(pos))
    marker_name = filename(1:4);
else
    marker_name = filename(pos(end)+1:pos(end)+4);
end
mn_len = min(60,length(marker_name));

if (flag_P1)
    code_id = 'P1';
else
    code_id = 'C1';
end

% date(:,1) = four_digit_year(date(:,1));

%write header
fprintf(fid_obs,'     2.10           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n');
fprintf(fid_obs,'goGPS                                                       PGM / RUN BY / DATE \n');
fprintf(fid_obs,'%-60sMARKER NAME         \n',marker_name(1:mn_len));
fprintf(fid_obs,'                                                            OBSERVER / AGENCY   \n');
fprintf(fid_obs,'                    %-20s                    REC # / TYPE / VERS \n', receiver);
fprintf(fid_obs,'                    %-20s                    ANT # / TYPE        \n', antenna);
fprintf(fid_obs,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', pos_R(1), pos_R(2), pos_R(3));
fprintf(fid_obs,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
fprintf(fid_obs,'     1     1                                                WAVELENGTH FACT L1/2\n');
fprintf(fid_obs,['     8    ' code_id '    L1    S1    D1    P2    L2    S2    D2      # / TYPES OF OBSERV \n']);
fprintf(fid_obs,'%10.3f                                                  INTERVAL            \n', interval);
fprintf(fid_obs,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
    date(1,1), date(1,2), date(1,3), date(1,4), date(1,5), date(1,6));
fprintf(fid_obs,'                                                            END OF HEADER       \n');

%-------------------------------------------------------------------------------

date(:,1) = two_digit_year(date(:,1));

%number of records
N = length(time_R);

%write data
for i = 1 : N
    
    sat = find(pr1_R(:,i) ~= 0);
    n = length(sat);
    
    %if no observations are available, do not write anything
    if (n > 0)
        fprintf(fid_obs,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
            date(i,1), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6), n);
        if (n>12)
            for j = 1 : 12
                fprintf(fid_obs,'G%02d',sat(j));
            end
            fprintf(fid_obs,'\n');
            fprintf(fid_obs,'%32s','');
            for j = 13 : n
                fprintf(fid_obs,'G%02d',sat(j));
            end
        else
            for j = 1 : n
                fprintf(fid_obs,'G%02d',sat(j));
            end
        end
        fprintf(fid_obs,'\n');
        for j = 1 : n
            %L1
            if (abs(pr1_R(sat(j),i)) > 1e-100)
                fprintf(fid_obs,'%14.3f %1d',pr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            else
                fprintf(fid_obs,'                ');
            end
            if (abs(ph1_R(sat(j),i)) > 1e-100)
                fprintf(fid_obs,'%14.3f %1d',ph1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            else
                fprintf(fid_obs,'                ');
            end
            fprintf(fid_obs,'%14.3f %1d',snr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            fprintf(fid_obs,'%14.3f %1d',dop1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            %L2
            if (abs(pr2_R(sat(j),i)) > 1e-100)
                fprintf(fid_obs,'%14.3f %1d',pr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
            else
                fprintf(fid_obs,'                ');
            end
            fprintf(fid_obs,'\n'); %end of line (5 observations per line)
            if (abs(ph2_R(sat(j),i)) > 1e-100)
                fprintf(fid_obs,'%14.3f %1d',ph2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
            else
                fprintf(fid_obs,'                ');
            end
            fprintf(fid_obs,'%14.3f %1d',snr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
            fprintf(fid_obs,'%14.3f %1d',dop2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
            fprintf(fid_obs,'\n');
        end
    end
end

%close RINEX observation file
fclose(fid_obs);

fprintf(['done.\n']);
