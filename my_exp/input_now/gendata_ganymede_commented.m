% =========================================================================
% gendata_ganymede.m
% -------------------------------------------------------------------------
% Builds the input/forcing files and edits the namelists ("data" and
% "data.shelfice") that MITgcm reads at runtime, for a simulation of
% Ganymede's subsurface ocean.
%
% MATLAB orientation notes (for newcomers):
%   - "%"  starts a single-line comment.
%   - "%%" starts a "section" that you can run with Ctrl+Enter.
%   - ";"  at the end of a line suppresses console output.
%   - "."  before *, /, ^ means element-wise (e.g.  A.*B  is element-wise
%          multiplication, while A*B is matrix multiplication).
%   - Array indexing uses 1-based parentheses: A(1) is the first element.
%   - "a:b:c" is a range from a to c in steps of b (MATLAB's "colon").
%   - sprintf builds a string; fprintf prints to console / a file.
%   - fopen/fwrite/fclose write raw binary files (used here so MITgcm can
%     read them directly).
%
% MITgcm orientation notes:
%   - The simulation reads two text namelists: "data" (core model) and
%     "data.shelfice" (the ice-shelf package).
%   - func_replace_string(filename, key, newline) finds a line containing
%     'key' in 'filename' and replaces it with 'newline'. That is how this
%     script edits the namelists in place.
%   - The binary inputs use big-endian ('b') single-precision ('real*4').
% =========================================================================


%% load packages -----------------------------------------------------------
% Add MITgcm's MATLAB utilities (e.g. rdmds for reading model output) and
% the GSW (Gibbs SeaWater) toolbox (used here to compute density, freezing
% temperature, thermal/haline expansion coefficients, etc.) to the path.
MITROOT='/home/cgliddon/'; % <-- where MITgcm is installed on your machine
addpath([MITROOT,'/MITgcm/utils/matlab/'])
addpath(genpath([MITROOT,'/MITgcm/my_exp/gsw/'])); % gsw must live here


%% initial-condition mode --------------------------------------------------
% restartTS controls how the initial T (potential temperature) and S
% (salinity) fields are produced.  0 = build analytical IC profiles from
% scratch.  Other modes restart from a previous run; see header comment.
restartTS=0;
restartTSpath='';     % directory containing previous run output (modes >0)
restartname='';       % suffix used when naming the IC binary files
                      % (used when restartTS ~= 0; see file-naming logic
                      % inside the "if restartTS<4" block below)
gridpath='./';        % grid files for pickup-restart modes (restartTS>3)
restartTSiter=0;      % which model iteration to restart from (NaN=time mean)


%% ice-surface temperature -------------------------------------------------
% Ts0 is the temperature at the very top of the ice shell (ice/space
% interface).  meridionalTs picks a latitudinal profile for it.
% Reference values: Ganymede ~120 K, Callisto ~126 K (Squyres 1980),
% Titan ~94 K (Jennings et al. 2016, nearly latitude-independent).
Ts0=120;            % [K]
meridionalTs=1;     % 0=constant, 1=Ojakangas & Stevenson 1989, 3=Nadeau & McGehee 2018
obl=0;              % obliquity [deg]; only matters for some Ts schemes


%% interior diffusivity and boundary mixing --------------------------------
% kappav is the *vertical* tracer diffusivity [m^2/s].  meridionalkappav
% lets it vary spatially; here it is uniform (=0).
meridionalkappav=0;  % 0=uniform, 1=equator-pole gradient, 2=different in BL
kappav_interior=0.001;
kappav_bnd=1;
kappav_eq=0.01;
kappav_pole=0.01;
kappav_width=30/2;   % [deg] e-folding scale of equator-pole gradient

% Transfer coefficients at the ice/ocean and core/ocean boundaries.
% Heat & salt coefficients have units of velocity [m/s] -- they convert a
% temperature/salinity contrast into a flux.  "drag" is a linear bottom-
% drag coefficient (also m/s).
iceshelf_heatcoef=1e-5;
iceshelf_saltcoef=1e-5;
iceshelf_dragcoef=1e-4;
bottom_heatcoef=1e-5;
bottom_saltcoef=1e-5;
bottom_dragcoef=1e-4;


%% domain ------------------------------------------------------------------
% Htot = full water-column thickness (ice + ocean) being modeled here.
% Note: Ganymede's ice/ocean partitioning is poorly constrained (Vance
% et al. 2014); see comment in source for context.  See Hice0 below for
% the ice thickness used here.
Htot=430e3;     % [m] total depth
nx=20;          % # longitude cells
ny=84;          % # latitude  cells
nr=86;          % # vertical (radial) levels
dx=2;           % [deg] zonal grid spacing
dy=2;           % [deg] meridional grid spacing
yyM=84;         % domain spans -yyM ... +yyM in latitude (deg)

% Bathymetry style: 0 = flat seafloor at depth Htot;
%                   1 = flat + a cosine^2(lat) bump (deltaH = 30 km).
bathyType=0;


%% equation of state -------------------------------------------------------
Sref=60;          % [g/kg] reference (mean) salinity, ~6% by weight
eos='MDJWF';      % nonlinear EOS (McDougall, Jackett, Wright, Feistel)
usePT=1;          % 1 => convert in-situ T to potential T inside shelfice code
% NOTE: per source TODO, the exact effect of usePT inside MITgcm's
% shelfice_thermodynamics.F should be confirmed.


%% ice-shell geometry ------------------------------------------------------
% Hice = Hice0 + Hice_P1*P1(lat) + Hice_P2*P2(lat) + Hice_P3*P3(lat),
% where P1,P2,P3 are Legendre polynomials of sin(lat).  Setting all P's
% to 0 (default below) gives a uniform-thickness shell.
Hice0=70e3;       % [m] mean shell thickness (Vance et al. 2018)
Hice_P1=0e3;      % asymmetric (hemispheric) component
Hice_P2=0e3;      % equator-pole component
Hice_P3=0;
realtopo=0;       % 0=use analytical Hice; 1=use zonal mean of file; 2=use 2D field
realtopopath='/home/wanying/Hemingway_Mittal_2019_Enceladus_nominal_shell_thickness_Fig11d/Enceladus_nominal_shell_thickness_Fig11d.tab';
Mice_randic=0*10.0;  % [m] amplitude of random ice-mass perturbation
                     % (set to literal "0*10.0" so you can flip the leading
                     % 0 to 1 to enable; equivalent to 0 right now)


%% heating partition and tidal heating shape -------------------------------
% Total heat loss through the ice (=Hcond, computed below) must be balanced
% by tidal heating in the ice and/or geothermal heating from the core.
% Htide0_portion = fraction supplied by ice-shell tidal heating; the rest
% comes from the bottom (core).
Htide0_portion=0;   % 0 = all heating from the core
twodtide=0;         % 1 = also include longitudinal (zonal) tidal pattern
qbotvary=1;         % 1 = use Beuthe (2019) poleward-amplified bottom heating
qbot0=-1;           % <0 => choose qbot0 to close the heat budget;
                    %  >0 => force qbot0 to that value [W/m^2]

% Spherical-harmonic coefficients of the tidal-dissipation pattern in the
% shell (relative to the Y00 mode), from Beuthe (2019) for Enceladus
% parameters.  Order is [Y20 Y40 Y22 Y42 Y44].
Htidemode=[0.250, 0.0825, -0.0834, -0.0546, -0.0562];
% Same idea for the "mix+bend" deformation mode; first entry here is the
% Y00 amplitude itself.
Hmixbendmode=[0.124,0.196,-0.0199,-0.0656,0.0132,0.0136];
addmixbend=1;       % include the mix+bend mode?

% Tidal dissipation and conductive heat loss scale as Hice^p:
ptide=-2.0;
pcond=-1.0;

% "Virtual" Legendre topography that affects only tidal/conductive terms
% (lets you decouple the heating geometry from the actual ice geometry).
HtidePs_ice=[0,0,0,0,0,0];
useHtidePsinHcond=0;
ptide_ext=-2.0;
pcond_ext=-1.0;

tiltheat=0.0;       % hemispheric asymmetry knob


%% ice freezing/melting evolution ------------------------------------------
PrescribeFreezing=0; % 1=use thin-shell ice flow model (Kang & Flierl 2020)
EvolveIce=0;         % 0=fix ice geometry, 1=let it evolve
SHIdtFactor=1;       % artificially accelerates ice evolution


%% physical constants ------------------------------------------------------
rhoice=917.;        % [kg/m^3] ice density
kappa0=651;         % [W/m] coefficient in T-dependent ice conductivity:
                    %       k(T) = kappa0/T  (so units are W/m, not W/m/K)
a0=2634.1e3;        % [m] Ganymede's mean radius
G=6.67e-11;         % [N m^2 / kg^2] Newton's constant
M=1.4819e23;        % [kg] Ganymede's mass

rhoout=1000;        % [kg/m^3] outer (ocean) reference density
% Solid-core density inferred so total mass equals M with a uniform-density
% ocean of thickness Htot above:
rhocore=3*M/(4*pi*(a0-Htot)^3) - rhoout*((a0/(a0-Htot))^3-1);

g0=G*M/a0^2;        % [m/s^2] surface gravity
eta_melt=3e14;      % [Pa s] ice viscosity at the melting point
eta_max_Q=Inf;      % cap on Q-dependent viscosity (none)
Ea=59.4e3;          % [J/mol] activation energy for ice creep
Rg=8.31;            % [J/mol/K] universal gas constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER-EDITABLE BLOCK ENDS HERE %%%%%%%%%%%%
% Below this point, derived quantities are computed and the namelists
% "data" and "data.shelfice" are edited via func_replace_string.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -- save a backup copy of the namelists ---------------------------------
% Useful in case this script clobbers them (running it is destructive).
system('cp -f data data_back');
system('cp -f data.shelfice data.shelfice_back');

% A short string used to name the binary files for non-default topography.
% NOTE the (slightly counter-intuitive) convention: "flat" => suffix
% present; non-flat => suffix is empty.
if (Hice_P1==0 && Hice_P2==0 && Hice_P3==0 && realtopo==0)
    appendix='flat';
else
    appendix='';
end


%% -- domain coordinates and Legendre polynomials -------------------------
% MITgcm needs cell-centre and cell-edge coordinates.  Here we build a
% regular lat-lon grid spanning [-yyM, yyM] in latitude.
nxy = nx*ny;
d2r = pi/180;

yc = -yyM+dy/2 : dy : yyM;             % cell-centre latitudes  (length ny)
xc = ((dx/2):dx:(nx*dx))' - nx*dx/2;   % cell-centre longitudes (length nx)
yi = -yyM:dy:yyM;                      % cell-edge   latitudes  (length ny+1)
xi = (0:dx:nx*dx)' - nx*dx/2;          % cell-edge   longitudes (length nx+1)

clat = cosd(yc);                       % cos(latitude)
slat = sind(yc);                       % sin(latitude)

% Vertical layer thicknesses [m].  Currently 86 uniform 5-km layers,
% summing to 430 km = Htot.  Two commented examples above show how you
% would build a graded grid (thinner near the ice, thicker below).
dh = 5.0e3*ones(1,86);
hf = [0, cumsum(dh)];                  % depths of layer interfaces (top->bot)

% Legendre polynomials P_n(sin(lat)) -- handy for projecting fields onto
% standard spherical-harmonic-like meridional modes.
P1 = slat;
P2 = (3/2).*slat.^2 - (1/2);
P3 = (5/2).*slat.^3 - (3/2).*slat;
P4 = (35*slat.^4 - 30.*slat.^2 + 3)/8.0;
P5 = (63.*slat.^5 - 70.*slat.^3 + 15.*slat)./8;
P6 = (231.*slat.^6 - 315.*slat.^4 + 105.*slat.^2 - 5)/16.0;
Ps = [P1; P2; P3; P4; P5; P6];

% Area weights for meridional averaging.  On a sphere, area scales like
% cos(lat); the form below uses the difference of sin(lat) at cell edges
% (which is the exact integral of cos), normalised to mean=1.
wgt1 = ny * (sind(yi(2:end))-sind(yi(1:end-1))) ./ (sind(yi(end))-sind(yi(1)));
wgt  = repmat(wgt1, [nx,1]);

prec = 'real*4';   % MITgcm reads single-precision big-endian floats


%% -- equation of state ----------------------------------------------------
% Compute reference density and EOS coefficients at the ice/ocean
% interface, where T equals the local freezing temperature.
P_inter0 = rhoice*g0*Hice0;                          % [Pa] pressure at top of ocean
Tfreeze0 = gsw_t_freezing(Sref, P_inter0/1e4);       % [degC] gsw uses dbar = 1e4 Pa
rhoNil   = gsw_rho   (Sref, Tfreeze0, P_inter0/1e4); % reference density
sBeta    = gsw_beta  (Sref, Tfreeze0, P_inter0/1e4); % haline contraction coeff
tAlpha   = gsw_alpha (Sref, Tfreeze0, P_inter0/1e4); % thermal expansion coeff
pKappa   = gsw_kappa_CT_exact(Sref, 0, P_inter0/1e4); % isentropic compress.
                                                      % (computed but not
                                                      % currently written
                                                      % to the namelist)

fprintf('\nrhoNil=%.5g\n', rhoNil)
fprintf('\nsBeta=%.4g, tAlpha=%.4g, pKappa=%.4g\n', sBeta, tAlpha, pKappa)

% Patch the values into the "data" namelist.  These values are still
% needed even with a nonlinear EOS because MITgcm uses rhoNil/tAlpha/sBeta
% as reference values for, e.g., the linearization of buoyancy in the
% momentum equations.
func_replace_string('data','rhoNil',  sprintf('rhoNil=%.3f,',  rhoNil))
func_replace_string('data','tAlpha',  sprintf('tAlpha=%.3g,',  tAlpha))
func_replace_string('data','sBeta',   sprintf('sBeta=%.3g,',   sBeta))
func_replace_string('data','eosType', sprintf('eosType=%s,',   eos))


%% -- gravity profile g(r) -------------------------------------------------
% Bulk density of the moon (using surface gravity & radius).  Inside a
% sphere of *uniform* density, g(r) is linear in r; here we account for a
% lighter outer shell of density rhoout and a denser core, which is why
% g grows toward the core for thick oceans.
rhobulk = g0/(G*4*pi/3*a0);
h = (hf(2:end)+hf(1:end-1))./2;        % depth at layer mid-points
g = g0*a0^2*(1 - rhoout/rhobulk*(a0^3-(a0-h).^3)./(a0^3))./(a0-h).^2;

func_replace_string('data','gravity=', sprintf('gravity=%f,',g0))
fprintf('gravity=%f', g0)
func_replace_string('data','rSphere=', sprintf('rSphere=%f,',a0))
fprintf('rSphere=%f\n', a0)

% Write g(r) to a binary file MITgcm can read.
fname = 'gravity_r_ganymede.bin';
fid   = fopen(fname,'w','b'); fwrite(fid,g,prec); fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data','gravityFile',['gravityFile=''',fname,''','])
fprintf(['gravityFile=''',fname,''','])


%% -- bathymetry ----------------------------------------------------------
% hh holds (negative) seafloor depth at each (x,y) cell.  Walls of zero
% depth at the north & south boundaries close off the domain.
hh           = ones(nx,ny);
hh(:,1)      = zeros(nx,1);
hh(:,ny)     = zeros(nx,1);
hh           = -Htot * hh;

if (bathyType==1)
    % Optional cosine^2(lat) bump in the seafloor (shallower at the equator).
    deltaH = 30e3;
    phis = linspace(-pi/2, pi/2, ny);
    h    = deltaH * cos(phis).^2;
    hh   = hh + repmat(h, nx, 1);
end

if (bathyType==0)
    file_name = ['flat_',sprintf('%d',Htot/1e3),'km.bin'];
else
    file_name = ['myBathy_',sprintf('%d',Htot/1e3),'km.bin'];
end
fid = fopen(file_name,'w','b'); fwrite(fid,hh,prec); fclose(fid);
fprintf(['\nwrite file: ', file_name])
func_replace_string('data','bathyFile',['bathyFile=''',file_name,''','])
fprintf(['bathyFile=''',file_name,''','])


%% -- ice topography (thickness Hice as a function of lat,lon) ------------
if (realtopo==0)
    % Analytical: mean + Legendre perturbations.
    Hice = Hice0 + Hice_P1.*P1 + Hice_P2.*P2 + Hice_P3.*P3;
    Hice = repmat(reshape(Hice,[1,ny]),[nx,1]);
else
    % Read 720 x 1440 (lat x lon) ice-thickness map from a binary file.
    fid   = fopen(realtopopath);
    Hice_ = fread(fid,'double','ieee-be');
    Hice_ = reshape(Hice_,[720,1440])*1e3;   % convert km -> m
    if realtopo==1
        % zonally average, then interpolate onto our latitude grid
        y_    = [-90:0.25:90-0.25]+0.25/2; y_ = [-90,y_,90];
        Hice_ = mean(Hice_,2);
        Hice_ = [Hice_(1);Hice_;Hice_(end)];
        Hice  = interp1(y_,Hice_,yc,'linear');
    else
        % full 2D bilinear interpolation onto our (xc,yc) grid
        y_    = [-90:0.25:90-0.25]+0.25/2; y_ = [-90,y_,90];
        x_    = [0:0.25:360];
        Hice_ = [Hice_, Hice_(:,1)];                  % wrap in longitude
        Hice_ = [Hice_(1,:); Hice_; Hice_(end,:)];    % pad in latitude
        [x2_,y2_] = meshgrid(x_,y_);
        [xc2,yc2] = meshgrid(xc,yc);
        Hice  = interp2(x2_,y2_,Hice_,xc2,yc2);
    end
    Hice0 = mean(mean(Hice.*wgt));     % overwrite Hice0 with new mean
end


%% -- ice-shelf "draft" file and ice-mass file ----------------------------
% Hunder = depth (m, negative) of the ice/ocean interface.  Treating the
% ice as exactly co-located with the surface gives Hunder = -Hice; the
% commented line above gives the proper Archimedean-floating depth and is
% the right thing to use if your model interprets Hunder that way.
Hunder = -Hice;

if isempty(appendix), fname='icetopo.bin';
else                  fname=['icetopo_',appendix,'.bin']; end
fid = fopen(fname,'w','b'); fwrite(fid,Hunder,prec); fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data.shelfice','SHELFICEtopoFile', ...
                    ['SHELFICEtopoFile=''',fname,''','])
fprintf(['SHELFICEtopoFile=''',fname,''','])

% Mass per unit area of the ice column.  Used by MITgcm to compute the
% pressure that the ice exerts on the ocean.  Optional random perturbation.
Miceshelf = Hice.*rhoice;
if Mice_randic
    Miceshelf_ = Miceshelf + rhoice.*smooth(rand(nx,ny)).*Mice_randic;
else
    Miceshelf_ = Miceshelf;
end

if isempty(appendix), fname='iceShelf_Mass.bin';
else                  fname=['iceShelf_Mass_',appendix,'.bin']; end
fid = fopen(fname,'w','b'); fwrite(fid,Miceshelf_,prec); fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data.shelfice','SHELFICEmassFile', ...
                    ['SHELFICEmassFile=''',fname,''','])
fprintf(['SHELFICEmassFile=''',fname,''','])


%% -- freezing temperature under spatially varying ice topography ---------
% MITgcm internally uses a *linear* freezing-T law of the form
%    Tf = fa0*S + fc0 + fb*P
% with S in g/kg and P in dbar.  Below we use the same coefficients so
% that the Tfreeze field we hand to the model is consistent with what the
% model itself computes.
fa0 = -0.0575;
fb  = -7.61e-4;
fc0 =  0.0901;

P_thin   = rhoice*g0*min(Hice(:));
P_thick  = rhoice*g0*max(Hice(:));
P_bot    = rhoout*g0*Htot;

% Freezing T from the full GSW (nonlinear) reference -- for diagnostics.
Tfreeze_thin  = gsw_t_freezing(Sref, P_thin/1e4);
Tfreeze_thick = gsw_t_freezing(Sref, P_thick/1e4);

% Freezing T using the model's linear law:
Tfreeze0_modelcode      = fa0*(Sref) + fc0 + fb*P_inter0/1e4;
Tfreeze_thin_modelcode  = fa0*(Sref) + fc0 + fb*P_thin /1e4;
Tfreeze_thick_modelcode = fa0*(Sref) + fc0 + fb*P_thick/1e4;

fprintf('\nTfreeze0=%.4g, Tfreeze_thin=%.4g,  Tfreeze_thick=%.4g', ...
        Tfreeze0, Tfreeze_thin, Tfreeze_thick)
fprintf('\nTfreeze0_mod=%.4g, Tfreeze_thin_mod=%.4g,  Tfreeze_thick_mod=%.4g', ...
        Tfreeze0_modelcode, Tfreeze_thin_modelcode, Tfreeze_thick_modelcode)

[Tfreeze0_SA, Tfreeze0_P] = gsw_t_freezing_first_derivatives(Sref, P_inter0/1e4);

% Reference 1D vertical T and S profiles for the "data" namelist
% (used by MITgcm to linearise about a hydrostatic background state).
% "nr*<value>" syntax means "repeat <value> nr times".
func_replace_string('data','tRef', sprintf('tRef=%d*%.3f,', nr, Tfreeze0_modelcode))
func_replace_string('data','sRef', sprintf('sRef=%d*%.3f,', nr, Sref))

% Optional pressure correction (currently set to 1 i.e. off).  The
% commented expression accounts for the fact that the pressure at the
% ice/ocean interface depends on radius via the gravity profile.
%p_correction=((1.-rhoout/rhobulk)./(1.-Hice/a0)+(rhoout/rhobulk).*(1.-Hice./a0/2.));
p_correction = 1;
P_inter      = Miceshelf .* g0 .* p_correction;
Tfreeze      = fa0*(Sref) + fc0 + fb*P_inter/1e4;   % 2D field (nx x ny)

if usePT
    % Convert in-situ Tfreeze to potential temperature (referenced to the
    % surface).  fPT_T is a helper function that should live somewhere on
    % the MATLAB path.
    Thfreeze = fPT_T(Tfreeze, P_inter/1e4, Sref);
    func_replace_string('data.shelfice','usePT','usePT=.TRUE.,')
    fprintf('usePT=.TRUE.')
else
    Thfreeze = Tfreeze;
    func_replace_string('data.shelfice','usePT','usePT=.FALSE.,')
    fprintf('usePT=.FALSE.')
end


%% -- conductive heat loss through the ice shell --------------------------
% First build the surface temperature Ts(lat).
Ts = Ts0.*ones(1,ny);
if meridionalTs==1
    % Ojakangas & Stevenson (1989): low latitudes follow the equilibrium
    % profile cos(lat)^(1/4); near the poles, with finite obliquity, an
    % asymptotic expansion is used.
    ycr     = yc.*d2r;
    ipolar  = abs(yc) > 90 - obl;
    Ts(~ipolar) = Ts0 .* (clat(~ipolar)).^(1/4);
    Ts( ipolar) = Ts0 .* (((pi/2-abs(ycr(ipolar))).^2 + (obl*d2r)^2)./2).^(1/8);
    % recalculate mean surface temperature
    Ts0=sum(Ts.*clat)/sum(clat);
elseif meridionalTs==3
    % Nadeau & McGehee (2018) Legendre expansion of insolation.
    flux_prof = 1 ...
        - (5/8)    .*legendreP(2,cosd(obl)).*legendreP(2,slat) ...
        - (9/64)   .*legendreP(4,cosd(obl)).*legendreP(4,slat) ...
        - (65/1024).*legendreP(6,cosd(obl)).*legendreP(6,slat);
    Ts_prof = flux_prof.^(1/4);
    Ts      = Ts_prof * Ts0;
end
func_replace_string('data.shelfice','meridionalTs', ...
                    sprintf('meridionalTs=%d,',meridionalTs))

Tssmt = reshape(smooth(Ts),[1,ny]);    % light smoothing for plotting
figure(1)
hold on
yyaxis right
plot(yc,Tssmt,'k-','LineWidth',2)
set(gca,'YColor','k','FontSize',16)
ylabel('Surface Temperature (K)')
xlabel('lat (deg)')

% Conductive heat flux through the shell.  Using k(T)=kappa0/T and
% integrating the steady 1D heat equation across the shell gives the
% closed-form Hcond = kappa0/Hice * ln(T_base/T_surface).
% Then a power-law thickness scaling Hice^pcond is applied.
Hcond = kappa0 .* log((Tfreeze0+273.15)./Ts) ./ Hice0;
Hcond = Hcond .* (Hice./Hice0).^pcond;

func_replace_string('data.shelfice','pcond=',     ['pcond=',    sprintf('%f',pcond),    ','])
fprintf('pcond=%f,',pcond)
func_replace_string('data.shelfice','pcond_ext=', ['pcond_ext=',sprintf('%f',pcond_ext),','])
fprintf('pcond_ext=%f,',pcond_ext)

if useHtidePsinHcond
    func_replace_string('data.shelfice','useHtidePsinHcond','useHtidePsinHcond=.TRUE.,')
    condtopo = 1;
    for i=1:length(HtidePs_ice)
        condtopo = condtopo + HtidePs_ice(i)*Ps(i,:);
        func_replace_string('data.shelfice', sprintf('HP%d_ice',i), ...
                            sprintf('HP%d_ice=%f,',i,HtidePs_ice(i)))
        fprintf('HP%d_ice=%f,',i,HtidePs_ice(i))
    end
    Hcond = Hcond .* (condtopo).^pcond_ext;
else
    func_replace_string('data.shelfice','useHtidePsinHcond','useHtidePsinHcond=.FALSE.,')
end

figure(1); hold on; yyaxis left
plot(yc,1e3.*Hcond,'g-.','LineWidth',2)
Hcond0 = mean(mean(Hcond.*wgt),2);     % global-mean conductive flux [W/m^2]
plot(yc,yc.*0+1e3.*Hcond0,'g-.')
fprintf('\nHcond=%.3g W/m^2', Hcond0)

func_replace_string('data.shelfice','SHELFICEkappa', ...
                    sprintf('SHELFICEkappa=%.1f,',-kappa0))
fprintf('SHELFICEkappa=%.1f,',-kappa0)
func_replace_string('data.shelfice','SHELFICEthetaSurface', ...
                    sprintf('SHELFICEthetaSurface=%.2f,',Ts0-273.15))
fprintf('SHELFICEthetaSurface=%.2f,',Ts0-273.15)
func_replace_string('data.shelfice','obliquity', ...
                    sprintf('obliquity=%.2f,',obl))
fprintf('obliquity=%.2f,',obl)


%% -- tidal heating in the ice shell --------------------------------------
% Total mean tidal heating that the model needs is whatever fraction
% Htide0_portion of the conductive loss is supposed to come from the
% shell (the rest will be supplied at the bottom).
Htidemean = Htide0_portion * Hcond0;

% Real surface spherical harmonics (un-normalised in convenient form).
one  = ones(nx,ny);
Y00  = sqrt(1/4/pi).*one;
Y20  = sqrt(5/4/pi)*(1.5*slat.^2-0.5).*one;
Y40  = sqrt(9/4/pi)*(35/8*slat.^4 - 30/8*slat.^2 + 3/8).*one;
if nx>1 && twodtide
    c2lon = cosd(2.*xc);
    c4lon = cosd(4.*xc);
    Y22 = sqrt(5/4/pi/24)   .* (3-3*slat.^2)                       .* (2*c2lon);
    Y42 = sqrt(9/4/pi/360)  .* (7.5.*(7*slat.^2-1).*(1-slat.^2))   .* (2*c2lon);
    Y44 = sqrt(9/4/pi/40320).* (105.*clat.^4)                      .* (2*c4lon);
else
    Y22 = Y00.*0;  Y42 = Y00.*0;  Y44 = Y00.*0;
end

% Membrane-mode dissipation pattern.  The sqrt(4*pi) renormalises so the
% global mean of the shape is 1 when only the Y00 term is present.
Htideprof = sqrt(4*pi).*( Y00 + Htidemode(1).*Y20 + Htidemode(2).*Y40 + ...
                          Htidemode(3).*Y22 + Htidemode(4).*Y42 + Htidemode(5).*Y44 );

% Optional virtual-topography modulation (HtidePs_ice).
tidetopo = 1;
for i=1:length(HtidePs_ice)
    tidetopo = tidetopo + HtidePs_ice(i)*Ps(i,:);
    func_replace_string('data.shelfice', sprintf('HP%d_ice',i), ...
                        sprintf('HP%d_ice=%f,',i,HtidePs_ice(i)))
    fprintf('HP%d_ice=%f,',i,HtidePs_ice(i))
end
Htideprof = Htideprof .* tidetopo.^ptide_ext;

if addmixbend
    func_replace_string('data.shelfice','addmixbend','addmixbend=.TRUE.,')
    fprintf('addmixbend=.TRUE.')
    Hmixbendprof = sqrt(4*pi).*( Hmixbendmode(1)*Y00 + Hmixbendmode(2).*Y20 + ...
                                 Hmixbendmode(3).*Y40 + Hmixbendmode(4).*Y22 + ...
                                 Hmixbendmode(5).*Y42 + Hmixbendmode(6).*Y44 );
else
    func_replace_string('data.shelfice','addmixbend','addmixbend=.FALSE.,')
    fprintf('addmixbend=.FALSE.')
    Hmixbendmode = 0;
    % Define an explicit zero mix+bend profile so it can still be added
    % into the Htidenorm expression below without erroring out.
    Hmixbendprof = zeros(size(Y00));
end

% Normalisation factor so that, after all the shape and thickness
% modulation, the area-mean tidal heat flux equals Htidemean.
Htidenorm = mean(mean( ((Hice./Hice0).^ptide.*Htideprof + Hmixbendprof) .* wgt ), 2);
Htide0    = Htidemean./Htidenorm;
fprintf('\nHtidemean=%.3g, Htide0=%.3g', Htidemean, Htide0)

func_replace_string('data.shelfice','Htide0', ['Htide0=',sprintf('%f',Htide0),','])
fprintf('Htide0=%f,',Htide0)
if EvolveIce
    func_replace_string('data.shelfice','Htide0_portion', ...
                        ['Htide0_portion=',sprintf('%f',Htide0_portion),','])
    fprintf('Htide0_portion=%f,',Htide0_portion)
else
    % If ice is fixed, comment out the line in the namelist (leading "#").
    func_replace_string('data.shelfice','Htide0_portion','#Htide0_portion=-1,')
    fprintf('Htide0_portion=%f,',Htide0_portion)
end

func_replace_string('data.shelfice','ptide=',    ['ptide=',    sprintf('%f',ptide),    ','])
fprintf('ptide=%f,',ptide)
func_replace_string('data.shelfice','ptide_ext=',['ptide_ext=',sprintf('%f',ptide_ext),','])
fprintf('ptide_ext=%f,',ptide_ext)

if twodtide
    func_replace_string('data.shelfice','tide2d','tide2d=.TRUE.,');  fprintf('tide2d=.TRUE.')
else
    func_replace_string('data.shelfice','tide2d','tide2d=.FALSE.,'); fprintf('tide2d=.FALSE.')
end

func_replace_string('data.shelfice','tiltheat', ...
                    ['tiltheat=',sprintf('%f',tiltheat),','])

figure(1); hold on; yyaxis left
plot(yc, Htideprof.*1e3.*Hcond0./Htidenorm, 'r-', 'LineWidth',2)
set(gca,'YColor','k','FontSize',16); box on
yl = ylim; ylim([0 yl(2)])
ylabel('Heat flux (mW/m^2)')


%% -- bottom (geothermal) heat flux ---------------------------------------
% If qbot0<0, choose qbot0 so the global heat budget closes:
%   surface conductive loss = ice tidal heating + bottom heating
% The (a0-...)^2/(a0-Htot)^2 factor accounts for the area difference
% between the top of the ocean (just below the ice) and the seafloor.
if qbot0<0
    qbot_inter = (Hcond0 - Htidemean);   % flux at TOP of ocean
    qbot0 = qbot_inter*(a0 - Hice0*(rhoice/rhoNil))^2 / (a0-Htot)^2;
end
qbot = qbot0*ones(nx,ny);

% Optional poleward-amplified profile from Beuthe (2019).  Normalised so
% area-mean = 1 (so changing the profile only redistributes heat).
if qbotvary
    qprofile = 1.08449 + 0.252257*cosd(2*(90-yc)) + 0.00599489*cosd(4*(90-yc));
    qprofile = qprofile/mean(qprofile.*wgt1);
    qbot     = qbot.*qprofile;
else
    qprofile = 1;
end
qbot = qbot(:);

file_name = sprintf('Q_bottom_%dmW.bin', round(qbot0*1e3));
fid = fopen(file_name,'w','b'); fwrite(fid,qbot,prec); fclose(fid);
fprintf(['\nwrite file: ', file_name])
func_replace_string('data','geothermalFile', ...
                    ['geothermalFile=''',file_name,''','])
fprintf(['geothermalFile=''',file_name,''','])

figure(1); hold on; yyaxis left
plot(yc,1e3.*qprofile.*Hcond0,'Color',[148,55,255]./256, ...
     'LineStyle','--','LineWidth',2)
xlim([yc(1), yc(end)])
saveas(gcf,'heat_profile.png')


%% -- surface and bottom masks (currently disabled) -----------------------
% These would be used by the RBCS package to nudge T or S in specific
% layers.  Not used in this configuration.


%% -- interface exchange rates and latent heat ----------------------------
Cp = 4000;   % [J/kg/K] ocean specific heat capacity (assumed constant)

% NOTE: this overwrites the 2D Tfreeze field with a Tfreeze "boosted" by
% the temperature contrast required to drive qbot_inter through the ice/
% ocean interface (Q = rho*Cp*gamma*deltaT, where gamma=iceshelf_heatcoef).
% The variable is reassigned but does not appear to be used downstream of
% here -- if you intended to use the boosted field somewhere, double-check.
Tfreeze = Tfreeze + qbot_inter/(iceshelf_heatcoef*rhoNil*Cp);

func_replace_string('data','HeatCapacity_Cp',  sprintf('HeatCapacity_Cp=%.3g,',Cp))
fprintf('\nHeatCapacity_Cp=%.3g,',Cp)
func_replace_string('data.shelfice','SHELFICEheatTransCoeff', ...
    sprintf('SHELFICEheatTransCoeff=%.3g,',iceshelf_heatcoef))
fprintf('SHELFICEheatTransCoeff=%.3g,',iceshelf_heatcoef)
func_replace_string('data.shelfice','SHELFICEsaltTransCoeff', ...
    sprintf('SHELFICEsaltTransCoeff=%.3g,',iceshelf_saltcoef))
fprintf('SHELFICEsaltTransCoeff=%.3g,',iceshelf_saltcoef)
func_replace_string('data','bottomDragLinear', ...
    sprintf('bottomDragLinear=%.3g,',bottom_dragcoef))
fprintf('bottomDragLinear=%.3g,',bottom_dragcoef)
func_replace_string('data.shelfice','SHELFICEDragLinear', ...
    sprintf('SHELFICEDragLinear=%.3g,',iceshelf_dragcoef))
fprintf('SHELFICEDragLinear=%.3g,',iceshelf_dragcoef)

SHELFICElatentHeat = 334000;   % [J/kg] latent heat of fusion of water ice
func_replace_string('data.shelfice','SHELFICElatentHeat', ...
    ['SHELFICElatentHeat=',sprintf('%f',SHELFICElatentHeat),','])


%% -- meridional / boundary-layer variation of vertical diffusivity -------
if meridionalkappav==1
    % Equator-pole gradient: kappa_v decays with |lat| with e-folding
    % scale kappav_width.
    kappav1 = kappav_pole + (kappav_eq-kappav_pole).*exp(-abs(yc)/kappav_width);
    kappav2 = repmat(kappav1,[nx,1]);
    kappav3 = repmat(kappav2(:),[1,nr]);
    fname   = sprintf('kappav3d_eq%.1gpole%.1g.bin',kappav_eq,kappav_pole);
    fprintf('kappa_v varying with latitude. kappav_pole=%f, kappav_eq=%f. writing file %s.\n', ...
            kappav_pole,kappav_eq,fname)
elseif meridionalkappav==2
    % Enhanced mixing inside a thin boundary layer near the ice/ocean
    % interface (decays exponentially with relative depth).
    Hunder1 = reshape(Hunder,[nxy,1]);
    relh    = max((Hunder1+h)./(Hunder1+Htot),0);
    kappav3 = kappav_interior + (kappav_bnd-kappav_interior).*exp(-relh.*nr./1);
    fname   = sprintf('kappav3d_bnd%.1ginterior%.1g.bin',kappav_bnd,kappav_interior);
    fprintf('kappa_v varying with latitude. kappav_interior=%f, kappav_bnd=%f. writing file %s.\n', ...
            kappav_interior,kappav_bnd,fname)
end
if meridionalkappav
    fid = fopen(fname,'w','b'); fwrite(fid,kappav3,prec); fclose(fid);
    func_replace_string('data','diffKrFile', ['diffKrFile=''',fname,''','])
    fprintf(['diffKrFile=''',fname,''','])
end


%% -- initial T and S fields ----------------------------------------------
% Goal: a stably stratified IC that is consistent with the ice/ocean
% boundary conditions (T at ice base = local freezing point) and that
% won't crash the model on step 1.

if restartTS>=0
if restartTS==0
    %---------------------------------------------------------------
    % Mode 0: build analytical IC from scratch.
    %---------------------------------------------------------------

    % Start with T = freezing T (column-by-column) at every depth.
    Tini1   = repmat(reshape(Thfreeze,[nxy,1]),[1,nr]);
    Hunder1 = reshape(Hunder,[nxy,1]);

    % Mean ice/ocean interface depth, treating the ice as Archimedean.
    Hunder0 = -(rhoice/rhoNil).*Hice0;

    % Locate the layer where the ice base sits, plus one extra layer for
    % the boundary-layer scheme used in shelfice.
    ibnd  = find(hf > -Hunder0, 1, 'first') + 1;
    h_bnd = hf(ibnd) + Hunder0 - 1.5*dh(ibnd-1)/2;

    % Relative depth (0 at top of BL, 1 at seafloor) and a slightly
    % "stretched" version that accounts for spherical geometry.
    relh  = max((Hunder1 + h - h_bnd)./(Hunder1 + Htot), 0);
    drelh = relh(:,2:end) - relh(:,1:end-1);
    drelh = drelh .* ((a0-mean(hf(2:end-1)))./(a0-hf(2:end-1))).^2.2;
    relh1 = [zeros(nxy,1), cumsum(drelh,2)];
    relh0 = max((h+Hunder0)./(Htot+Hunder0),0);
    relh  = 0*relh0 + 1*relh1;        % effectively relh = relh1
                                      % (written this way so you can blend)

    Tini = Tini1;

    % Pick the sign of the imposed top-to-bottom T contrast.  At high
    % salinity (>15) cold water is denser, so we want bottom warmer than
    % top (delTemp0 < 0 in our convention).  Below 15 g/kg the EOS is
    % such that warm-fresh water is denser, so the sign flips.
    if Sref>15
        delTemp0   = min(Thfreeze(:)) - max(Thfreeze(:));
        delTempvar = delTemp0;
    else
        delTemp0   = max(Thfreeze(:)) - min(Thfreeze(:));
        delTempvar = -delTemp0;
    end

    delTemp_pattern = delTemp0 - delTempvar*cosd(2.*yc);
    delTemp_pattern = repmat(delTemp_pattern,[nx,1]);
    delTemp_pattern = delTemp_pattern(:);
    delTemp = delTemp_pattern;       % comment in source: ~10 mK by trial

    Tini = Tini + delTemp.*ones(nxy,nr).*relh;

    dTemp = round(delTemp0*1000);    % in milli-K, used in filename
    sfx   = sprintf('%i%s',dTemp,'mK.bin');

    % Add small symmetric noise mid-column (kicks off convection).
    dTnoise = 0.01;
    dTn = rand([nx,ny]);
    dTn(:,ny/2+1:ny) = dTn(:,ny/2:-1:1);   % mirror to keep N-S symmetric
    dTn = dTn - mean(mean(dTn,1).*wgt1);   % zero mean
    dTn = dTn(:);
    Tini(:,nr/2) = Tini(:,nr/2) + dTnoise*dTn;

    if isempty(appendix), fTname=['theta_ini.',sfx];
    else                  fTname=['theta_ini_',appendix,'.',sfx]; end

    % --- salinity: uniform Sref everywhere (delSalt0=0 => no perturbation)
    Sini     = Sref + zeros(nxy,nr);
    delSalt0 = 0;
    Sini = Sini + min(h-abs(Hunder0),0)/(abs(Hunder0)-min(abs(Hunder(:))))*delSalt0;
    sfx  = sprintf('%i%s',delSalt0*1e3,'mpsu.bin');

    if isempty(appendix), fSname=['salt_ini.',sfx];
    else                  fSname=['salt_ini_',appendix,'.',sfx]; end

elseif restartTS<4
    %---------------------------------------------------------------
    % Modes 1-3: restart from the binary T/S/grid output of a previous run.
    %---------------------------------------------------------------
    [Tini,its,~] = rdmds([restartTSpath,'/T'], restartTSiter);
    [Sini,its,~] = rdmds([restartTSpath,'/S'], restartTSiter);
    xcrest = squeeze(mean(rdmds([restartTSpath,'/XC'],NaN),1));
    ycrest = squeeze(mean(rdmds([restartTSpath,'/YC'],NaN),1));
    rcrest = squeeze(mean(rdmds([restartTSpath,'/RC'],NaN),1));

    % If the user asked for a time-mean (NaN in restartTSiter), rdmds
    % returns several time slices; collapse them.
    if length(its)>1
        Tini = mean(Tini,4); Tini = Tini(:,:,:,1);
        Sini = mean(Sini,4); Sini = Sini(:,:,:,1);
    end

    if restartTS==1
        % horizontal-mean profile -> tile to (nx,ny)
        Tini = mean(Tini,[1,2]); Tini = repmat(Tini,[nx,ny,1]);
        Sini = mean(Sini,[1,2]); Sini = repmat(Sini,[nx,ny,1]);
    end

    if restartTS==2
        % zonal mean (-> 2D field), then interpolate onto current grid.
        disp('entering 2')
        Tini = squeeze(mean(Tini,1));
        Sini = squeeze(mean(Sini,1));

        if size(Tini,1)~=ny || size(Tini,2)~=nr
            % step 1: extrapolate into bathymetry-zeroed cells by walking
            % differences upward from the bottom (avoids NaN in interp2)
            disp('fill in')
            Tini_ = Tini; Tini_(Tini_==0) = NaN;
            Tbase = Tini_(:,end);
            Tdiff = -diff(Tini_,1,2); Tdiff(isnan(Tdiff))=0;
            Tini(:,1:end-1) = Tbase + cumsum(Tdiff,2,'reverse');
            Tini(1,:)   = Tini(2,:);
            Tini(end,:) = Tini(end-1,:);

            Sini_ = Sini; Sini_(Sini_==0) = NaN;
            Sbase = Sini_(:,end);
            Sdiff = -diff(Sini_,1,2); Sdiff(isnan(Sdiff))=0;
            Sini(:,1:end-1) = Sbase + cumsum(Sdiff,2,'reverse');
            Sini(1,:)   = Sini(2,:);
            Sini(end,:) = Sini(end-1,:);

            % step 2: spline interp onto our (lat, depth) grid
            disp('interpolate')
            [rrrest,yyrest] = meshgrid(rcrest,ycrest);
            [rr,yy]         = meshgrid(-h,yc);
            Tini_ = Tini; Tini = interp2(rrrest,yyrest,Tini_,rr,yy,'spline');
            Sini_ = Sini; Sini = interp2(rrrest,yyrest,Sini_,rr,yy,'spline');
        end

        Tini = repmat(reshape(Tini,[1,ny,nr]),[nx,1,1]);
        Sini = repmat(reshape(Sini,[1,ny,nr]),[nx,1,1]);
    end

elseif restartTS==4
    %---------------------------------------------------------------
    % Mode 4: take an MITgcm netCDF "pickup" file from another run and
    % use it as the *initial* condition (iteration 0) for this one.
    % expandpickups is a helper that re-tiles the data onto your grid.
    %---------------------------------------------------------------
    if ~exist('./pickup.0000000000.t001.nc','file')
        system('rm -rf initial_condition')
        system('mkdir initial_condition')
        system(['cp ',gridpath,'*.nc ./initial_condition/'])
        system(['rename 0000000001 ',sprintf('%010d',0),' ./initial_condition/*.nc'])
        expandpickups(restartTSpath,restartTSiter,'./initial_condition/',0)
        system('cp -f initial_condition/pickup*.nc .')
    end
    func_replace_string('data','pickupSuff','pickupSuff=''0000000000'',')

elseif restartTS==5
    %---------------------------------------------------------------
    % Mode 5: continue another run from iteration restartTSiter
    % (the pickup is renamed to keep the iteration counter).
    %---------------------------------------------------------------
    if ~exist(['./pickup.',sprintf('%010d',restartTSiter),'.t001.nc'],'file')
        system('rm -rf restart_condition')
        system('mkdir restart_condition')
        system(['cp ',gridpath,'/*.nc ./restart_condition/'])
        system(['rename 0000000000 ',sprintf('%010d',restartTSiter),' ./restart_condition/*.nc'])
        expandpickups(restartTSpath,restartTSiter,'./restart_condition/',restartTSiter)
        system('cp -f restart_condition/pickup*.nc .')
    end
    func_replace_string('run.sub','export currentiter', ...
                        ['export currentiter=',sprintf('%d',restartTSiter)])

elseif restartTS==6
    %---------------------------------------------------------------
    % Mode 6: convergence trick.  Use the difference between two earlier
    % runs (long-term tendency) to nudge the current pickup forward.
    % (Likely to be removed in a future iteration of this script.)
    %---------------------------------------------------------------
    disp('start interpolate nc file')
    restartTSiterout = restartTSiter - 1;
    if ~exist(sprintf('./pickup.%010d.t001.nc',restartTSiterout),'file')
        system('rm -rf restart_condition')
        system('mkdir restart_condition')
        system(['cp ',gridpath,'*.nc ./restart_condition/'])
        system(['rename ',sprintf('%010d',restartTSiter),' ',sprintf('%010d',restartTSiterout),' ./restart_condition/*.nc'])

        T_0 = mean(mean(rdmds(['run',sprintf('%d',diffiter0),'/T'],NaN),4),1);
        S_0 = mean(mean(rdmds(['run',sprintf('%d',diffiter0),'/S'],NaN),4),1);
        T_1 = mean(mean(rdmds(['run',sprintf('%d',diffiter1),'/T'],NaN),4),1);
        S_1 = mean(mean(rdmds(['run',sprintf('%d',diffiter1),'/S'],NaN),4),1);

        dTIC = repmat(reshape((T_1-T_0),[1,ny,nr]),[nx,1,1]);
        dSIC = repmat(reshape((S_1-S_0),[1,ny,nr]),[nx,1,1]);
        expandpickups(restartTSpath,restartTSiter,'./restart_condition/', ...
                      restartTSiterout,1,dTIC,dSIC)
        system('cp -f restart_condition/pickup*.nc .')
    end
end % end of restartTS dispatch

if restartTS<4
    % Add a tiny amount of noise to the bottom layer (kick-starts
    % bottom-driven convection without violating any constraints).
    Tini = reshape(Tini,[nxy,nr]);
    Sini = reshape(Sini,[nxy,nr]);
    dTnoise = 0.001e-6;
    dTn = rand([nxy,nr]);
    Tini(:,nr) = Tini(:,nr) + dTnoise.*dTn(:,nr);

    fTname = ['Tini_restart_',restartname,'.bin'];
    fSname = ['Sini_restart_',restartname,'.bin'];

    fid = fopen(fTname,'w','b'); fwrite(fid,Tini,prec); fclose(fid);
    fprintf(['\nwrite file: ',fTname]);
    func_replace_string('data','hydrogThetaFile', ...
                        ['hydrogThetaFile=''',fTname,''','])
    fprintf(['hydrogThetaFile=''',fTname,''','])

    fid = fopen(fSname,'w','b'); fwrite(fid,Sini,prec); fclose(fid);
    fprintf(['\nwrite file: ',fSname]);
    func_replace_string('data','hydrogSaltFile', ...
                        ['hydrogSaltFile=''',fSname,''','])
    fprintf(['hydrogSaltFile=''',fSname,''','])
end
end % if restartTS>=0


%% -- ice flow (creep) coefficient ----------------------------------------
% Computes a depth-integrated, viscosity-averaged "flow constant" for a
% Newtonian, Arrhenius-temperature-dependent ice rheology, used by the
% shelfice package to translate slope into ice flux.  The integral2 call
% does a 2D numerical integration over (T', T) where T' < T (a triangle in
% the temperature plane); the integrand is 1/eta(T') * ln(Tm/T') / T' / T.
Tsmean = mean(Ts.*wgt1);
Tm     = Tfreeze0 + 273.15;       % melting T at the ice base [K]

ice_flow0 = 2*((rhoNil-rhoice)*rhoice/rhoNil*g0)*Hice0^3 ...
            /((log(Tm/Tsmean))^3) ...
            * integral2( @(Tp,T) ...
                  1./min(eta_melt.*exp(-Ea/Rg/Tm + (Ea/Rg)./Tp), eta_max_Q) ...
                  .* ((log(Tm./Tp)))./Tp./T.*(Tp<=T), ...
                  Tsmean, Tm, Tsmean, Tm );

SHI_iceflow = ice_flow0 / Hice0^3 / rhoice^3;

func_replace_string('data.shelfice','SHI_iceflow', ...
                    ['SHI_iceflow=',sprintf('%g',SHI_iceflow),','])
fprintf('SHI_iceflow=%g,',SHI_iceflow)
func_replace_string('data.shelfice','SHELFICElatentHeat', ...
                    ['SHELFICElatentHeat=',sprintf('%g',SHELFICElatentHeat),','])


%% -- ice-evolution flags -------------------------------------------------
% These swap the model between "rigid lid" (ice geometry frozen) and
% "evolving ice" (mass-stepping + remeshing turned on).  Note the side
% effects: enabling evolving ice forces quasi-hydrostatic + nonlinear
% free-surface and turns OFF non-hydrostatic mode.
func_replace_string('data.shelfice','PrescribeFreezing', ...
                    sprintf('PrescribeFreezing=%d,',PrescribeFreezing))
fprintf('PrescribeFreezing=%d\n',PrescribeFreezing)

if EvolveIce
    func_replace_string('data','quasiHydrostatic','quasiHydrostatic=.TRUE.,')
    func_replace_string('data','nonlinFreeSurf','nonlinFreeSurf=4,')
    func_replace_string('data','nonHydrostatic','nonHydrostatic=.FALSE.,')
    func_replace_string('data.shelfice','SHIdtFactor', ...
        ['SHIdtFactor=',sprintf('%f',SHIdtFactor),','])
    func_replace_string('data.shelfice','SHELFICEMassStepping', ...
        'SHELFICEMassStepping=.TRUE.,')
    func_replace_string('data.shelfice','SHELFICERemeshFrequency', ...
        'SHELFICERemeshFrequency=2592000.,')   % every 30 days
    func_replace_string('data.shelfice','SHELFICEwriteState', ...
        'SHELFICEwriteState=.TRUE.,')
    fprintf('SHIdtFactor=%f\n',SHIdtFactor)
else
    func_replace_string('data','quasiHydrostatic','quasiHydrostatic=.FALSE.,')
    func_replace_string('data','nonlinFreeSurf','nonlinFreeSurf=0,')
    func_replace_string('data','nonHydrostatic','nonHydrostatic=.TRUE.,')
    func_replace_string('data.shelfice','SHELFICEMassStepping', ...
        'SHELFICEMassStepping=.FALSE.,')
    func_replace_string('data.shelfice','SHELFICERemeshFrequency', ...
        'SHELFICERemeshFrequency=0,')
    func_replace_string('data.shelfice','SHELFICEwriteState', ...
        'SHELFICEwriteState=.TRUE.,')
end