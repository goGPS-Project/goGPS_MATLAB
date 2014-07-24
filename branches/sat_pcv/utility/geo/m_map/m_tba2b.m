function m_tba2b(fnam);
% M_TBA2B Converts the ASCII TerrainBase 5-minute bathymetry database
%         (size 56Mb) available from
%           ftp://ncardata.ucar.edu/datasets/ds759.2/tbase.Z
%         into a binary file of 2-byte integers that can be read by
%         M_TBASE to provide high-resolution global bathymetry.
%
%         To use this file, first
%
%         a) get and uncompress the tbase.Z file from the above URL into the
%            current directory.
%
%         b) run this function:
%
%            m_tba2b(PATHNAME)
%
%         to store the resulting binary (of size 18Mb) as PATHNAME/tbase.int
%
%         c) Edit the PATHNAME setting in M_TBASE to point to the
%            location of this file.
%
%         d) delete the ASCII file tbase.
%

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Oct/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

if nargin==0,
 fnam='.';
end;

fnam=[fnam '/tbase.int'];

fid=fopen('tbase','rt');

if fid==-1,
 error('Cannot find file called ''tbase'' ');
end;

fidb=fopen(fnam,'w');
if fidb==-1,
 error(['Cannot open file ''' fnam '''']);
end;

for k=1:466560,
 data=fscanf(fid,'%6d',20);
 fwrite(fidb,data,'int16');
 if rem(k,2000)==0,
   disp([ int2str(k) '/466450 lines processed']);
 end;
end;


