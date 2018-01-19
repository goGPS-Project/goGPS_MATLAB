% M_Map - mapping toolbox (Author: rich@eos.ubc.ca)
% Version 1.4i  Nov 2017
%
% You have collected your data, loaded it into Matlab, analyzed 
% everything to death, and now you want to make a simple map showing 
% how it relates to the world. 
%
% But you can't. 
%
% Instead you have to figure out how to save all your data, and 
% then read it into a mapping program, and then spend all that extra 
% time figuring out why the mapping program doesn't give you what
% you expected it would...
%
% No more! 
%
%                            Announcing M_Map v1.4! 
%
% M_Map is a set of mapping tools written for Matlab v5. These include: 
%
%    1. Routines to project data in 18 different spherical 
%       projections (and determine inverse mappings) 
%    2. A grid generation routine to make nice axes with 
%       limits either in long/lat terms or in planar
%       X/Y terms. 
%    3. A coastline database (with 1/4 degree resolution) 
%    4. A global elevation database (1 degree resolution) 
%    5. Hooks into freely available high-resolution coastlines and
%       bathymetry/topography.
%
% M_Map v1.4 is available via the web at 
%
%       http://www.eos.ubc.ca/~rich/
%
%
% Toolbox contents
%
%    Contents.m    - This file
%    m_demo.m      - demonstrates a few different maps.
%
%  User-callable functions
%
%    m_proj.m      - initializes projections
%    m_coord.m     - converts between geomagnetic and geographic coords.
%
%    m_grid.m      - draws grids 
%    m_scale       - forces map to a given scale.
%    m_ruler       - draw a scale ruler
%
%    m_ungrid.m    - erases map elements (if you want to change parameters)
%
%    m_coast.m     - draws a coastline
%    m_elev.m      - draws elevation data from 1 degree database
%
%    m_tbase.m     - draws elevation data from 5-minute TerrainBase database
%    m_gshhs.m     - draws coastline from GSHHS with specified resolution
%    m_gshhs_c.m   - draws coastline from GSHHS crude database
%    m_gshhs_l.m   - draws coastline from GSHHS low-resolution database
%    m_gshhs_i.m   - draws coastline from GSHHS intermediate-resolution database
%    m_gshhs_h.m   - draws coastline from GSHHS high-resolution database
%    m_gshhs_f.m   - draws coastline from GSHHS full database
%    m_plotbndry.m - draws a political boundary from the DCW 
%    m_usercoast.m - draws a coastline using a user-specified subset database.
%
%    m_plot.m      - draws line data in map coords
%    m_line.m      - draws line data in map coords
%    m_text.m      - adds text data in map coords
%    m_legend.m    - draws a legend box
%    m_quiver.m    - draws arrows for vector data
%    m_contour.m   - draws contour lines for gridded data
%    m_contourf.m  - draws filled contours
%    m_patch.m     - draws patch data
%    m_pcolor.m    - draws pcolor data
%    m_streamline.m- draws streamlines
%    m_scatter.m   - draws scatter plot
%    m_annotation.m- annotation lines/boxes/text
%
%    m_track.m     - draws annotated tracklines
%    m_hatch.m     - hatched or speckled patches.
%    m_range_ring.m- draws range rings (spherical coords)
%    m_ellipse.m   - draws tidal ellipses (most requested ocean feature!)
%
%    m_ll2xy.m     - converts from long/lat to map coordinates
%    m_xy2ll.m     - converts from map coordinates to long/lat
%
%    m_geo2mag.m     - converts from long/lat to geomagnetic coordinates
%    m_mag2geo.m     - converts from geomagnetic coordinates to long/lat
%
%    m_lldist      - spherical distance/geodesics between points (long/lat coordinates)
%    m_xydist      - spherical distance between points (map projection coordinates)
%
%    m_fdist       - ellipsoidal geodesic forward calculation 
%    m_idist       - ellipsoidal geodesic inverse calculation 
%    m_geodesic    - points along ellipsoidal geodesics
%
%    m_tba2b.m     - used in installing high-resolution elevation database.
%
%    m_vec.m       - fancy arrows
%    m_windbarb.m  - barbed wind arrows
%
%    m_contfbar.m  - draws colorbars for contourf plots
%    m_colmap.m    - useful perceptually uniform colourmaps.
%    m_shaperead.m - reads ESRI shapefiles
%    mygrid_sand2.m- reads Sandwell and Smith bathymetry file
%
%    wysiwyg.m     - Sets figure window to match size/aspect of printed output
%
%  Internal functions (not meant to be user-callable)
%
%    private/mp_azim.m   - azimuthal projections
%    private/mp_cyl.m    - cylindrical projections (equatorial)
%    private/mp_conic.m  - conic projections
%    private/mp_tmerc.m  - transverse cylindrical projections
%    private/mp_utm.m    - elliptical universal transverse cylindrical projections
%    private/mp_omerc.m  - oblique cylindrical projection
%
%    private/mu_util.m   - various utility routines
%    private/mu_coast.m  - routines to handle coastlines.
%
%    private/mc_coords.m - coordinate systems based on different poles.
%    private/mc_ellips.m - parameters of different ellipsoidal earth models
%
%    private/m_coasts.mat- low-res coastline data
%
%  HTML documentation
%
%    map.html           - Home page, examples
%    private/mapug.html - User's guide
%    private/*gif       - examples.
%  
%
% Questions or problems; email me - rich@eos.ubc.ca.
%
% Rich Pawlowicz
% Dept. of Earth, Ocean, Atmospheric Sciences, Univ. of British Columbia, 
% 6339 Stores Rd., Vancouver, B.C. CANADA V6T 1Z4
% email: rich@eos.ubc.ca 
%
%
 

    

