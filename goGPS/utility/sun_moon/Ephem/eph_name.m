function [ fn ] = eph_name( denum, year )
%eph_name creates a naming convection for the file names of the binary ephemerides

  fn = Ephem.asc_name( denum, year );
  if ~isempty(fn)
    fn(1:3) = 'eph';
  end
end

