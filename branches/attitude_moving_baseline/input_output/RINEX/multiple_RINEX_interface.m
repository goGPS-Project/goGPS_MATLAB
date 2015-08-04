function [filename_obs] = multiple_RINEX_interface(filename_R_obs, filename_M_obs, mode)

%if multiple RINEX files in input
if (iscell(filename_R_obs))
    %if mode is not multi-receiver, fall back to the first rover filename
    if (~goGNSS.isMR(mode))
        filename_obs{1,1} = [filename_R_obs{1,1} filename_R_obs{1,2}];
        fprintf('... WARNING: multiple rover RINEX files in input for a single rover mode; using the first file.\n');
    else
        nFiles = length(filename_R_obs)-1;
        filename_obs = cell(nFiles,1);
        for f = 1 : nFiles
            filename_obs{f,1} = [filename_R_obs{1} filename_R_obs{f+1}];
        end
    end
else
    filename_obs{1} = filename_R_obs;
end

%if relative positioning
if (goGNSS.isDD(mode))
    filename_obs = [filename_obs; filename_M_obs];
end
