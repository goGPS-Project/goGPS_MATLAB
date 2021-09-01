function urlOut = urldecode(urlIn) 
% matlab encode/decode functions are not included in matlab runtime
% compiler (windows)

urlOut = char(java.net.URLDecoder.decode(urlIn,'UTF-8')); 

end 
