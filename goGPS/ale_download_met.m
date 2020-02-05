file_name = './met_list.lst';
out_dir = 'tmp_aria';

fid = fopen(file_name, 'wb');
for doy = 1 : 365
    remote_path = sprintf('ftp://igs.bkg.bund.de/IGS/obs/2019/%03d/m0se%03d0.19m.Z', doy, doy  );
    fprintf(fid, '%s\n', remote_path);
end
fclose(fid);
try
    mkdir(out_dir);
catch ex
end
dos(sprintf('/usr/local/bin/aria2c --input-file="%s" -d %s -j 8', file_name, out_dir))