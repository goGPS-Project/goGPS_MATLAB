function varargout = tracking( varargin )
%tracking  Track anonymized usage data
%
%  tracking(p,v,id) tracks usage to the property p for the product version
%  v and identifier id.  No personally identifiable information is tracked.
%
%  r = tracking(...) returns the server response r, for debugging purposes.
%
%  tracking('on') turns tracking on.  tracking('off') turns tracking off.
%  tracking('query') returns the tracking state.

%  tracking('spoof') sets the tracking settings -- domain, language,
%  client, MATLAB version, operating system version -- to spoof values.
%  tracking('reset') sets the tracking settings to normal values.
%
%  [t,s] = tracking('query') returns the tracking state t and settings s.

%  Copyright 2016 The MathWorks, Inc.
%  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $

persistent STATE USERNAME DOMAIN LANGUAGE CLIENT MATLAB OS
if isempty( STATE )
    STATE = getpref( 'Tracking', 'State', 'on' );
    if strcmp( STATE, 'snooze' ) % deprecated
        setpref( 'Tracking', 'State', 'on' )
        STATE = 'on';
    end
    if ispref( 'Tracking', 'Date' ) % deprecated
        rmpref( 'Tracking', 'Date' )
    end
    USERNAME = getenv( 'USERNAME' );
    reset()
end % initialize

switch nargin
    case 1
        switch varargin{1}
            case {'on','off'}
                STATE = varargin{1};
                setpref( 'Tracking', 'State', varargin{1} ) % persist
            case 'spoof'
                spoof()
            case 'reset'
                reset()
            case 'query'
                varargout{1} = STATE;
                varargout{2} = query();
            otherwise
                error( 'tracking:InvalidArgument', ...
                    'Valid options are ''on'', ''off'' and ''query''.' )
        end
    case 3
        switch nargout
            case 0
                if strcmp( STATE, 'off' ), return, end
                uri = 'https://www.google-analytics.com/collect';
                track( uri, varargin{:} );
            case 1
                uri = 'https://www.google-analytics.com/debug/collect';
                varargout{1} = track( uri, varargin{:} );
            otherwise
                nargoutchk( 0, 1 )
        end
    otherwise
        narginchk( 3, 3 )
end % switch

    function reset()
        %reset  Set normal settings
        
        DOMAIN = lower( getenv( 'USERDOMAIN' ) );
        LANGUAGE = char( java.util.Locale.getDefault() );
        CLIENT = getpref( 'Tracking', 'Client', uuid() );
        MATLAB = matlab();
        OS = os();
        
    end % reset

    function spoof()
        %spoof  Set spoof settings
        
        DOMAIN = randomDomain();
        LANGUAGE = randomLanguage();
        CLIENT = randomClient();
        MATLAB = randomMatlab();
        OS = randomOs();
        
    end % spoof

    function s = query()
        %query  Return settings
        
        s.Username = USERNAME;
        s.Domain = DOMAIN;
        s.Language = LANGUAGE;
        s.Client = CLIENT;
        s.Matlab = MATLAB;
        s.Os = OS;
        
    end % query

    function varargout = track( uri, p, v, s )
        %track  Do tracking
        
        a = sprintf( '%s/%s (%s)', MATLAB, v, OS );
        if isdeployed()
            ds = 'deployed';
        elseif strcmp( DOMAIN, 'mathworks' )
            ds = DOMAIN;
        else
            ds = 'unknown';
        end
        pv = {'v', '1', 'tid', p, 'ua', escape( a ), 'ul', LANGUAGE, ...
            'cid', CLIENT, 'ht', 'pageview', ...
            'dp', sprintf( '/%s', s ), 'ds', ds};
        [varargout{1:nargout}] = urlread( uri, 'Post', pv );
        
    end % track

end % tracking

function s = randomDomain()
%randomDomain  Random domain string

switch randi( 4 )
    case 1
        s = 'mathworks';
    otherwise
        s = hash( uuid() );
end

end % randomDomain

function s = randomLanguage()
%randomLanguage  Random language string

lo = java.util.Locale.getAvailableLocales();
s = char( lo(randi( numel( lo ) )) );

end % randomLanguage

function s = randomClient()
%randomClient  Random client identifier

s = uuid();

end % randomClient

function s = matlab()
%matlab  MATLAB version string

v = ver( 'MATLAB' );
s = v.Release;
s(s=='('|s==')') = [];

end % matlab

function s = randomMatlab()
%randomMatlab  Random MATLAB version string

releases = {'R2014b' 'R2015a' 'R2015b' 'R2016a' 'R2016b'};
s = releases{randi( numel( releases ) )};

end % randomMatlab

function s = os()
%os  Operating system string

if ispc()
    s = sprintf( 'Windows NT %s', ...
        char( java.lang.System.getProperty( 'os.version' ) ) );
elseif isunix()
    s = 'Linux x86_64';
elseif ismac()
    s = sprintf( 'Macintosh; Intel OS X %s', ...
        strrep( char( java.lang.System.getProperty( 'os.version' ) ), ' ', '_' ) );
else
    s = 'unknown';
end

end % os

function s = randomOs()
%randomOs  Random operating system string

switch randi( 3 )
    case 1
        versions = [5.1 5.2 6 6.1 6.2 6.3 10];
        s = sprintf( 'Windows NT %.1f', ...
            versions(randi( numel( versions ) )) );
    case 2
        s = 'Linux x86_64';
    case 3
        s = sprintf( 'Macintosh; Intel OS X 10_%d', ...
            randi( [10 12] ) );
end

end % randomOs

function s = escape( s )
%escape  Escape string

s = char( java.net.URLEncoder.encode( s, 'UTF-8' ) );

end % escape

function h = hash( s )
%hash  Hash string
%
%  See also: rptgen.hash

persistent MD5
if isempty( MD5 )
    MD5 = java.security.MessageDigest.getInstance( 'MD5' );
end

MD5.update( uint8( s(:) ) );
h = typecast( MD5.digest, 'uint8' );
h = dec2hex( h )';
h = lower( h(:) )';

end % hash

function s = uuid()
%uuid  Unique identifier

s = char( java.util.UUID.randomUUID() );

end % uuid