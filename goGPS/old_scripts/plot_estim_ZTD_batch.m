if (exist(sprintf('%s_ZTD.txt', file_name_base), 'file'))
    ZTD_batch = load(sprintf('%s_ZTD.txt', file_name_base));
    plot(reshape(ZTD_batch',1,size(ZTD_batch,1)*2880))
end