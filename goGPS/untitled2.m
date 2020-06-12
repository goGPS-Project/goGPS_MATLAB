% Copy file

source_dir = '/Volumes/Data/Volumes/WorkDisk/localArchive/${YYYY}/${DOY}/${I}';
dest_dir = '/Volumes/Data/ArchiveGNSS/STA/${YYYY}/${DOY}/${I}';
clc
for y = 2019:2020
    for d = 1:365
        for i = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            source_tmp = strrep(source_dir, '${YYYY}', sprintf('%04d', y));
            dest_tmp = strrep(dest_dir, '${YYYY}', sprintf('%04d', y));
            source_tmp = strrep(source_tmp, '${DOY}', sprintf('%03d', d));
            dest_tmp = strrep(dest_tmp, '${DOY}', sprintf('%03d', d));
            source_tmp = strrep(source_tmp, '${I}', i);
            dest_tmp = strrep(dest_tmp, '${I}', i);
            if exist(source_tmp, 'dir')
                fprintf('%4d %03d %c\n', y, d, i);
                if exist(dest_tmp, 'dir')
                    dos(['cp -f ' source_tmp, '/* ' dest_tmp '/']);
                else
                    dos(['cp -R ' source_tmp, ' ' dest_tmp '']);
                end
            end
        end
    end
end