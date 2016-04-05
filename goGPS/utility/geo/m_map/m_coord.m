function coordsys=m_coord(coordsys);
% M_COORD Initializes the coordinate system for varous conversions.
%
%   M_COORD('set') tells you the current coordinate system
%   M_COORD('get') gives you the current possibilities
%   M_COORD(SYSTEM) sets the coordinate system to SYSTEM.
%
%   Currently the available coordinate systems are:
%              'geographic' (usual lat/long)
%              'geomagnetic' (referenced to magnetic poles)

% Rich Pawlowicz (rich@ocgy.ubc.ca) nOV/2001
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


global MAP_COORDS

if nargin==0, coordsys='usage'; end;

coordsys=lower(coordsys);

coords=mc_coords('name');

switch coordsys,

   case 'get',
      disp(' ');
      disp('Available coordinate systems are:'); 
      for k=1:length(coords.name),
        disp(['     ' coords.name{k}]);
      end;
   
   case 'set',
      if isempty(MAP_COORDS),
         disp('No coordinate system initialized');
         m_coord('usage');
      else
         if nargout==0,
           disp(MAP_COORDS.name);
	 else
	   coordsys=MAP_COORDS;
	 end;    
      end;
 
   case 'usage',
      disp(' ');
      disp('Possible calling options are:');
      disp('  ''usage''                    - this list');
      disp('  ''set''                      - list of coordinate systems');
      disp('  ''get''                      - get current coordinate (if defined)');
      disp('  ''system''                   - initialize coordinate system\n');
   
   otherwise
     k=strmatch(coordsys,lower(coords.name));
     MAP_COORDS=mc_coords('parameters',coords.name{k});
   
end;   
    % Check to see if a non-lat/long coordinate system is being used.
