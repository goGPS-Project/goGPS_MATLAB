function urlOut = urlencode(urlIn) 
% matlab encode/decode functions are not included in matlab runtime
% compiler (windows)

urlOut = char(java.net.URLEncoder.encode(urlIn,'UTF-8'));

end 

