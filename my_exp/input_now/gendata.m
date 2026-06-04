% =========================================================================
% gendata.m
% -------------------------------------------------------------------------
% This function builds the the input/forcing files and edits the namelists 
% ("data", "data.shelfice", "data.hpimm") that MITgcm reads at runtime, 
% for a simulation of subsurface icy moon oceans.
% 
% To run (must use Matlab 2023a or later): matlab -batch "gendata"
%
% MATLAB orientation notes (for newcomers):
%   - "."  before *, /, ^ means element-wise (e.g.  A.*B  is element-wise
%          multiplication, while A*B is matrix multiplication).
%   - Array indexing uses 1-based parentheses: A(1) is the first element.
%   - "a:b:c" is a range from a to c in steps of b (MATLAB's "colon").
%   - sprintf builds a string; fprintf prints to console / a file.
%   - fopen/fwrite/fclose write raw binary files (used here so MITgcm can
%     read them directly).
%
% MITgcm orientation notes:
%   - func_replace_string(filename, key, newline) finds a line containing
%     'key' in 'filename' and replaces it with 'newline'. That is how this
%     script edits the namelists in place.
%   - The binary inputs use big-endian ('b') single-precision ('real*4').
%
% Code authors: Coleman Gliddon (cgliddon@mit.edu), Wanying Kang 
%   (wanying@mit.edu)
% =========================================================================

function gendata()

%% load packages
% Add MITgcm's MATLAB utilities (e.g. rdmds for reading model output) and
% the GSW (Gibbs SeaWater) toolbox (used here to compute density, freezing
% temperature, thermal/haline expansion coefficients, etc.) to the path.
MITROOT='/home/cgliddon/'; % <-- where MITgcm is installed on your machine
addpath([MITROOT,'/MITgcm/utils/matlab/'])
eos = eosMethods(); % defined in the file eosMethods.m, in same directory as gendata.m

%% set basic model parameters and configuration

% physical constants
G=6.674e-11;           % Newton's gravitational constant, N m^2 kg^-2
rho_Ih=917.;           % ice-Ih density, kg m^-3
k_Ih=651;              % ice-Ih thermal conductivity coefficient, W m^-1
kappa_wh=1.0e-3;       % water horizontal thermal diffusivity, m^2 s^-1
kappa_wr=1.0e-3;       % water vertical thermal diffusivity, m^2 s^-1
nu_wh=1.;              % water horizontal viscosity, m^2 s^-1
nu_wr=1.;              % water vertical viscosity, m^2 s^-1
alpha_T=4.5e-4;        % water thermal expansion coefficient, 
cpw=4.184e+03;         % water heat capacity, J kg^-1
pref=0.0;              % reference pressure for PT-T interconversion, Pa

% basic planetary parameters
satellite_name='myMoon';
M=1.4819e+23;          % planetary mass, kg
a=2634.1e+03;          % planetary radius, m
rotationPeriod=86400*7.15455; % planetary rotation period, s
obl=0.0;               % obliquity, degrees
Ts0=110;               % surface temperature, K
meridionalTs=0;        % 0: uniform surface temperature; 1: Ojakangas & Stevenson (1989) profile; 3: Nadeau & McGehee (2018) profile
                       % Note that Ts0 has different meanings depending on this choice: mean temp. for 0 & 3, equatorial temp. for 1
hice0=50.0e3;          % ice shell mean thickness
del_ht=1.0e3;          % topographic variation scale for top boundary, m
del_hb=1.0e3;          % topographic variation scale for bottom boundary, m
delv_theta = 0.0;      % superadiabatic ocean temperature gradient
S0=10.0;               % ocean salinity, g kg^-1

% field initialization flags
Tini_field=1;          % 0: isentropic; 1: superadiabatic following Gastine et al. scaling; ; 2: read from file
Sini_field=0;		   % 0: isohaline; 2: read from file
Vini_field=0;		   % 0: does nothing; 2: read from file

% model configuration flags
simpleTopo=true;       % true: topography assumed to vary as cos^2 of latitude; false: manually prescribed topography (overrides del_ht, del_hb)
addFloorice=false;     % false: no floor ice; true: high-pressure floor ice
bottom_ice_phase="V";  % reference bottom ice phase
% The following flags control whether tidal heating is included in the model
% and, if it is, how it's handled. 
% > useShellTides=0: no tidal heating
% > useShellTides=1: manually specified tidal heating
% > useShellTides=2: automatically computed tidal heating following Beuthe 
useShellTides=0;
useOceanTides=0;       % 0/1: excludes/includes ocean tidal heating (TODO: MUST BE SET TO ZERO FOR NOW)
Htide0_portion=0;      % fraction of heat flux into ice shell due to tidal heating (used only if useShellTides=2)
prescribeHcore=false;  % false: core heating assumed uniform with latitude; true: specify pattern manually
prescribe_delv_theta=true;  % false: delv_theta assumed fixed; true: delv_theta calculated iteratively from balance conditions

% model grid: number of cells, spacing
nx = 96; ny = 96;
dx = 1.75; dy = 1.75;
yyM = 84;                             % max. model latitude in both hemispheres
dz_max = 2.5e+03; dz_min = 8.0e+02;   % max and min grid spacing for vertical levels
drchfrac = 0.10;                      % how much do we allow the vertical spacing to change by for adjacent cells?

prec='real*4';                        % data save-out precision

%% coordinates, Legendre polynomials, and spherical harmonics
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
c2lon = cosd(2.*xc);
c4lon = cosd(4.*xc);

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

% Real surface spherical harmonics (un-normalised in convenient form).
one  = ones(nx,ny);
Y00  = sqrt(1/4/pi).*one;
Y20  = sqrt(5/4/pi)*(1.5*slat.^2-0.5).*one;
Y40  = sqrt(9/4/pi)*(35/8*slat.^4 - 30/8*slat.^2 + 3/8).*one;
Y22 = sqrt(5/4/pi/24)   .* (3-3*slat.^2)                       .* (2*c2lon);
Y42 = sqrt(9/4/pi/360)  .* (7.5.*(7*slat.^2-1).*(1-slat.^2))   .* (2*c2lon);
Y44 = sqrt(9/4/pi/40320).* (105.*clat.^4)                      .* (2*c4lon);

%% topography calculation

% set Legendre polynomial coefficients for top and bottom topographies
if simpleTopo
    hice_Pc  = [0.; -del_ht; 0.; 0.; 0.; 0.];
    bathy_Pc = [0.; -del_hb; 0.; 0.; 0.; 0.];
else
    % TODO: modify this section 
    hice_Pc  = [0.; 0.; 0.; 0.; 0.; 0.];
    bathy_Pc = [0.; 0.; 0.; 0.; 0.; 0.];  
end

% construct meridional profiles
hice  = hice0 + sum(hice_Pc.*Ps);
bathy = sum(bathy_Pc.*Ps);      % this is just anomaly from the mean ocean depth for now

%% freezing point and Hcond calculation

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
        - (5/8)    .*legendreP(2,cosd(obl)).*P2 ...
        - (9/64)   .*legendreP(4,cosd(obl)).*P4 ...
        - (65/1024).*legendreP(6,cosd(obl)).*P6;
    Ts_prof = flux_prof.^(1/4);
    Ts      = Ts_prof * Ts0;
end

% next, calculate mean freezing point at ice-ocean interface 
p1=(rho_Ih*G*hice0)/(a^2-a*hice0)*(M+2*pi*rho_Ih*a^2*hice0-2*pi*rho_Ih*a*hice0^2/3);
T1=eos.eval_freezing_point(p1,S0,"Ih");
theta1=eos.t_to_theta(p1,T1,S0,pref);    % potential temperature

% calculate conductive heat flux profile
pcond=-1.0;
Hcond = k_Ih .* log(T1./Ts) ./ hice0;
Hcond = Hcond .* (hice./hice0).^pcond;
Hcond0 = sum(Hcond.*clat)/sum(clat);    % mean heat flux

% TODO: include correction which incorporates tides

%% Htide calculation

if useShellTides==0
    Htide0 = 0;
elseif useShellTides>1
    % Spherical-harmonic coefficients of the tidal-dissipation pattern in the
    % shell (relative to the Y00 mode), from Beuthe (2019) for Enceladus
    % parameters.  Order is [Y20 Y40 Y22 Y42 Y44].
    % TODO: make this more general
    Htidemode=[0.250, 0.0825, -0.0834, -0.0546, -0.0562];
    % Same idea for the "mix+bend" deformation mode; first entry here is the
    % Y00 amplitude itself.
    Hmixbendmode=[0.124,0.196,-0.0199,-0.0656,0.0132,0.0136];
    addmixbend=1;       % include the mix+bend mode?

    if useShellTides==1
    % Total heat loss through the ice (=Hcond, computed below) must be balanced
    % by tidal heating in the ice and/or geothermal heating from the core.
    % Htide0_portion = fraction supplied by ice-shell tidal heating; the rest
    % comes from the bottom (core).
        Htidemean = Htide0_portion * Hcond0;
    else
        % TODO: insert your mean tidal heat flux here
        Htidemean = 0.0; 
    end

    ptide=-2.0;
    % Membrane-mode dissipation pattern.  The sqrt(4*pi) renormalises so the
    % global mean of the shape is 1 when only the Y00 term is present.
    Htideprof = sqrt(4*pi).*( Y00 + Htidemode(1).*Y20 + Htidemode(2).*Y40 + ...
                          Htidemode(3).*Y22 + Htidemode(4).*Y42 + Htidemode(5).*Y44 );

    if addmixbend
        Hmixbendprof = sqrt(4*pi).*( Hmixbendmode(1)*Y00 + Hmixbendmode(2).*Y20 + ...
                                    Hmixbendmode(3).*Y40 + Hmixbendmode(4).*Y22 + ...
                                    Hmixbendmode(5).*Y42 + Hmixbendmode(6).*Y44 );
    else
        Hmixbendmode = 0;
        % Define an explicit zero mix+bend profile so it can still be added
        % into the Htidenorm expression below without erroring out.
        Hmixbendprof = zeros(size(Y00));
    end

    % Normalisation factor so that, after all the shape and thickness
    % modulation, the area-mean tidal heat flux equals Htidemean.
    Htidenorm = mean(mean( ((hice./hice0).^ptide.*Htideprof + Hmixbendprof) .* wgt ), 2);
    Htide0    = Htidemean./Htidenorm;
end

%% Iterative calculation of self-consistent ocean depth
Hcore0_inter = Hcond0 - Htide0;
power = Hcore0_inter*4*pi*(a-hice0)^2;

% define function for nonlinear solving
function poc = calc_ocean_pres(r, rho_oc, G, a, hice0, M, rho_Ih, p1)
% Calculates the interior ocean pressure at radius r assuming uniform
% density rho_oc.
    r_inner = a - hice0;
    poc = rho_oc*G*(r_inner-r)/(r_inner*r) * (M - 4*pi/3*rho_Ih*(3*a*r_inner + hice0^2) ...
            - 2*pi/3*rho_oc*r_inner*(r_inner - r)*(2*r_inner + r));
    poc = poc + p1;
end

if Tini_field==0 || Tini_field==1
    % iteratively estimate rho_oc and H_oc
    H_iters = 5;
    rho_iters = 3;
    rho_oci = 1100.0;
    hoc_i   = 60.0e+03;
    trap_res = 31;

    for hi = 1:H_iters
        r_coords = linspace(a - hice0 - hoc_i, a - hice0, trap_res);  % Use hoc_i, not Hoc
        rho_p = @(p) eos.rho_pThetaS(p, theta1, S0, pref);
        
        for rhoi = 1:rho_iters
            % Pass parameters explicitly
            pres_r = @(r) calc_ocean_pres(r, rho_oci, G, a, hice0, M, rho_Ih, p1);
            integrand = zeros(size(r_coords));
            for idx = 1:length(r_coords)
                integrand(idx) = rho_p(pres_r(r_coords(idx)));
            end
            % Use hoc_i consistently
            rho_oci = trapz(r_coords, integrand) / hoc_i;
        end
        ac = a - hice0 - hoc_i;

        delta_theta = @(p) eos.t_to_theta(p, eos.eval_freezing_point(p, S0, bottom_ice_phase), S0, pref) - delv_theta - theta1;
        p_guess = calc_ocean_pres(a - hice0 - hoc_i, rho_oci, G, a, hice0, M, rho_Ih, p1);
        fprintf('Pressure guess %.6f yields freezing point:\n', p_guess);
        fprintf('%.6f\n', theta1 + delv_theta + delta_theta(p_guess));
        
        % Use single initial guess (fzero will use secant/quasi-Newton methods)
        p_freeze = fzero(delta_theta, p_guess);
        % fprintf('Freezing pressure: %.6f\n', p_freeze);
        
        % invert pres_func to find the freezing radius
        delta_pres = @(r) calc_ocean_pres(r, rho_oci, G, a, hice0, M, rho_Ih, p1) - p_freeze;
        r_guess1 = a - hice0;
        r_guess2 = 0.1 * a;
        hoc_i = a - fzero(delta_pres, [r_guess1, r_guess2]) - hice0;
        % fprintf('Ocean depth after %d iterations: %.6f\n', mi, hoc_i);
        if ~prescribe_delv_theta
            % recalculates delv_theta using Gastine et al. (2015) scaling
            % TODO: include multiple scalings for different rotation rates
	        delv_theta=sqrt(nu_wr*(ac/(a-hice0))*power/(alpha_T*rho_oci*cpw*g_top*hoc_i^2));
        end
        fprintf('Superadiabatic temperature gradient (K): %.6f\n',delv_theta)
    end
end
hoc0 = hoc_i; rho0 = rho_oci; p2=p_freeze;

% correct core heating mean value for change in surface area from ocean top to bottom
Hcore0 = Hcore0_inter*(a - hice0)^2/(a - hice0 - hoc0)^2;

%% prescribe bottom heat flux variations
% TODO: write up this section

if prescribeHcore
    fprintf('Warning: prescribeHcore not implemented yet \n')
else
    Hcore = repmat(Hcore0, [1, ny]);
end

%% create vertical grid
bathy_z = bathy - (hice0 + hoc0);     % bathymetry in z-coordinates
hliq = -min(hice) - min(bathy_z);     % total liquid depth of model

% calculate vertical grid spacing
ratio = 1.0 + drchfrac;
transition = [];
dz_current = dz_min;
while dz_current < dz_max
    transition = [transition; dz_current];
    dz_current = dz_current * ratio;
end
    
boundary_depth = 2.0 * sum(transition);
if boundary_depth > hliq
    error(sprintf( ...
        ['The two boundary transition layers together require ', ...
         '%.2f m, which exceeds H = %.2f m. ', ...
         'Decrease dz_min, increase dz_max, or increase H.'], ...
    boundary_depth, hliq));
end
    
interior_depth = hliq - boundary_depth;
n_interior = round(interior_depth / dz_max);

% Adjust dz_max to fill interior depth exactly
if n_interior > 0
    dz_interior = interior_depth / n_interior;
else
    % No interior cells — rescale transition layer uniformly to fill H
    scale = hliq / boundary_depth;
    dz_interior = [];
    transition = transition * scale;
end
    
% Assemble full column
if ~isempty(dz_interior)
    interior = repmat(dz_interior, n_interior', 1);
    dz = [transition; interior; flipud(transition)];
else
    dz = [transition; flipud(transition)];
end
dh = transpose([hice0; dz]);
fprintf("dh profile (km): %.6f\n", dh/1.0e+3);
hf = [0,cumsum(dh)];
hc=(hf(2:end)+hf(1:end-1))./2;  % centered depth coordinates
rc=a-hc;
nr=length(hc);

%% gravity profile

g0 = G*M/a^2;
g_top = p1/(rho_Ih*dh(1));
g = [g_top, 4*G*pi/3.*((rho_cr - rho0).*ac^3./rc(2:nr).^2 + rho0.*rc(2:nr))];

%% -- ice flow (creep) coefficient ----------------------------------------
% Computes a depth-integrated, viscosity-averaged "flow constant" for a
% Newtonian, Arrhenius-temperature-dependent ice rheology, used by the
% shelfice package to translate slope into ice flux.  The integral2 call
% does a 2D numerical integration over (T', T) where T' < T (a triangle in
% the temperature plane); the integrand is 1/eta(T') * ln(Tm/T') / T' / T.

eta_melt=3e14;      % [Pa s] ice viscosity at the melting point
eta_max_Q=Inf;      % cap on Q-dependent viscosity (none)
Ea=59.4e3;          % [J/mol] activation energy for ice creep
Rg=8.31;            % [J/mol/K] universal gas constant

ice_flow0 = 2*((rho0-rho_Ih)*rho_Ih/rho0*g0)*hice0^3 ...
            /((log(T1/Ts0))^3) ...
            * integral2( @(Tp,T) ...
                  1./min(eta_melt.*exp(-Ea/Rg/T1 + (Ea/Rg)./Tp), eta_max_Q) ...
                  .* ((log(T1./Tp)))./Tp./T.*(Tp<=T), ...
                  Ts0, T1, Ts0, T1 );

SHI_iceflow = ice_flow0 / Hice0^3 / rhoice^3;

% TODO: adapt this so the *full* ice flow calculation is in gendata.m 

%% initial dynamic-variable profiles
% temperature
if Tini_field == 0 || Tini_field == 1
    % linear temperature profile
    Tprof = theta1*ones(1, nr);
    Tprof(2:end)=Tprof(2:end)+delv_theta*(hc(2:end)-mean(hc(2:end)))/(hf(end)-hf(2));
    T_3d=repmat(Tprof, nxy, 1);
    % add random noise
    delTemp = 1.0e-3;
    dTn=2.0*rand(size(T_3d))-1.0;
    T_3d=T_3d + delTemp*dTn - 273.15;
end

% salinity (TODO)
% velocity (TODO)

%% write to namelists

eosType='SEAFRZ';
func_replace_string('data','rSphere=',sprintf('rSphere=%f,',a))
func_replace_string('data','eosType=',sprintf('eosType=%s,',eosType))
func_replace_string('data','rotationPeriod=',sprintf('rotationPeriod=%f',rotationPeriod))
func_replace_string('data','viscAh=',sprintf('viscAh=%f,',nu_wh))
func_replace_string('data','viscAr=',sprintf('viscAr=%f,',nu_wr))
func_replace_string('data','diffKhT=',sprintf('diffKhT=%f,',kappa_wh))
func_replace_string('data','diffKrT=',sprintf('diffKrT=%f,',kappa_wr))
func_replace_string('data','tAlpha=',sprintf('tAlpha=%f,',alpha_T))
func_replace_string('data','rhoConst=',sprintf('rhoConst=%f,',rho0));
func_replace_string('data','gravity=',sprintf('gravity=%f,',g0));

func_replace_string('data.shelfice','meridionalTs',sprintf('meridionalTs=%d,',meridionalTs))
func_replace_string('data.shelfice','SHI_iceflow', ['SHI_iceflow=',sprintf('%g',SHI_iceflow),','])

func_replace_string('data.hpimm', 'hpimm_addFloorice=', sprintf('hpimm_addFloorice=.%s.', addFloorice?"TRUE":"FALSE"))
func_replace_string('data.hpimm', 'hpimm_flooriceType=', sprintf('hpimm_flooriceType=%s', bottom_ice_phase))



%% write binary files

% gravity profile
fname=strcat('gravity_', satellite_name,'.bin');
fid=fopen(fname,'w','b');
fwrite(fid,g,prec);
fclose(fid);
func_replace_string('data','gravityFile=',['gravityFile=''',fname,''','])

% topography files
Miceshelf=hice.*rho_Ih;
fname='iceShelf_Mass.bin';
fid=fopen(fname,'w','b');
fwrite(fid,Miceshelf,prec);
fclose(fid);
func_replace_string('data.shelfice','SHELFICEmassFile', ['SHELFICEmassFile=''',fname,''','])

Hunder=repmat(-hice, [nx, 1]);
fname='icetopo.bin';
fid=fopen(fname,'w','b');
fwrite(fid,Hunder,prec);
fclose(fid);
func_replace_string('data.shelfice','SHELFICEtopoFile', ['SHELFICEtopoFile=''',fname,''','])

bathy_z = repmat(bathy_z, [nx, 1]);
fname=['my_bathy_',sprintf('%3d', -min(bathy_z)/1.0e+3),'km.bin'] ;
fid=fopen(fname,'w','b'); 
fwrite(fid,bathy_z,prec); 
fclose(fid);
func_replace_string('data','bathyFile=',['bathyFile=''',fname,''','])

% heat conduction
Hcond = repmat(Hcond, [nx,1]);
fname = strcat('Hcond_', satellite_name,'.bin');
fid=fopen(fname,'w','b');
fwrite(fid,Hcond,prec);
fclose(fid);

% tidal heating
Htide = repmat((hice./hice0).^ptide.*Htideprof + Hmixbendprof, [nx, 1]);
fname = strcat('Htide_', satellite_name,'.bin');
fid=fopen(fname,'w','b');
fwrite(fid,Htide,prec);
fclose(fid);

% bottom heating
Hcore = repmat(Hcore, [nx, 1]);
fname = strcat('Hcore_', satellite_name,'.bin');
fid=fopen(fname,'w','b');
fwrite(fid,Hcore,prec);
fclose(fid);
func_replace_string('data','geothermalFile=', ['geothermalFile=''',fname,''','])
% potential temperature
fname=strcat('theta_ini_', sprintf('%i%s',delTemp*1.0e3,'mK.bin'));
fid=fopen(fname,'w','b');
fwrite(fid,T_3d,prec);
fclose(fid);
func_replace_string('data','hydrogThetaFile',['hydrogThetaFile=''',fname,''','])



end