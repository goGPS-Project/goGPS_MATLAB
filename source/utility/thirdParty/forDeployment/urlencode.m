function urlOut = urlencode(urlIn)
%URLENCODE Replace special characters with escape characters URLs need

% Copyright 1984-2008 The MathWorks, Inc.

urlOut = char(java.net.URLEncoder.encode(urlIn,'UTF-8'));