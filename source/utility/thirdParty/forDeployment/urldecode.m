function urlOut = urldecode(urlIn)
%URLDECODE Replace URL-escaped strings with their original characters

% Copyright 1984-2008 The MathWorks, Inc.

urlOut = char(java.net.URLDecoder.decode(urlIn,'UTF-8'));