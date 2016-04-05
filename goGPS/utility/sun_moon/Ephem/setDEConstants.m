function [s,error] = setDEConstants( s, denum )
% sets constant values that do not depend on time
  if nargin < 2
    denum = s.de_number;
  end
  if isnumeric(denum)
    denum = sprintf('%d',denum);
  end
  s.de_number = 0;
  %   Set initializations and default values.

  %   Set the value of the record length according to what JPL ephemeris is
  %   being opened.
  [ st, error ] = Ephem.readHeader( denum );
  if error == 0
    s.RECORD_LENGTH = 8*st.NCOEFF;
    s.JPLAU = st.AU;
    s.EM_RATIO = st.EMRAT;
    s.masses(Ephem.Mercury) = st.GM1;
    s.masses(Ephem.Venus) = st.GM2;
    s.masses(Ephem.Mars) = st.GM4;
    s.masses(Ephem.Jupiter) = st.GM5;
    s.masses(Ephem.Saturn) = st.GM6;
    s.masses(Ephem.Uranus) = st.GM7;
    s.masses(Ephem.Neptune) = st.GM8;
    s.masses(Ephem.Pluto) = st.GM9;
    s.masses(Ephem.Sun) = st.GMS;
    s.masses(Ephem.EarthMoonBarycenter) = st.GMB;
    s.masses(Ephem.Moon) = st.GMB/(1+s.EM_RATIO);
    s.masses(Ephem.Earth) = s.masses(Ephem.EarthMoonBarycenter)-s.masses(Ephem.Moon);
    m = 0;
    for i=Ephem.Mercury:Ephem.Sun
      m = m + s.masses(i);
    end
    s.masses(Ephem.SolarSystemBarycenter) = m;
    s.IPT = st.IPT;
    s.de_number = denum;
  end
end
