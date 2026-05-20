% Updated version of Wanying Kang's gendata.m code for planetary structure and MITgcm initialization.
function gendata()
%% load packages
MITROOT='/home/cgliddon/'; % MITROOT should point to the directory where MITgcm is installed XX.
addpath([MITROOT,'/MITgcm/utils/matlab/'])
eos = eosMethods(); % defined in the file eosMethods.m, in same directory
%% Constants and free parameters
G = 6.674e-11;           % Newton's gravitational constant, N m^2 kg^-2
rho_Ih = 917.;           % ice-Ih density, kg m^-3
M = 1.0759e+23;          % planetary mass, kg
a = 2410.3e+03;          % planetary radius, m
Hice = 150.e+03;         % ice crust thickness, m
S = 10.0;                % salinity, g kg^-1
pref = 1.0e+05;          % reference pressure for PT-T interconversion
bottom_ice_phase = "V";  % reference bottom ice phase
dHbot = 30.0e3;          % topographic variation scale for bottom boundary
simpleBathy=1;           % 1: bathymetry is just -P2(sin(phi))
Tini_profile=2;          % 0: isentropic; 1: follows freezing point of varying bottom topography; 2: read from file
Sini_profile=2;		 % 0: isohaline; 2: read from file
init_Vel=1;		 % 0: does nothing; 1: read from file
prescribe_Qbot=0;	 % 0: does nothing; 1: generates and writes out a bottom heating pattern
satellite_name='myMoon';
restart_datfile='exp_data.mat';
restart_data = load(restart_datfile);

nx=336.;ny=336.;
dx=0.5;dy=0.5; yyM=84.;
dz_max = 2.0e+03; dz_min = 5.0e+02;   % max and min grid spacing for vertical levels
drchfrac = 0.10;                      % how much do we allow the vertical spacing to change by?

prec='real*4';                        % data save-out precision

%% start replacing parameters in data file
eosType='SEAFRZ';
func_replace_string('data','rSphere=',sprintf('rSphere=%f,',a))
func_replace_string('data','eosType',sprintf('eosType=%s,',eosType))



%% Obtain density estimates for 3-layer interior structure
% pressure at ocean surface
p0_ice = rho_Ih*G*Hice/(a - Hice)*(M/a + 2*pi*rho_Ih*Hice*(Hice/3 - a));
% freezing PT at ocean surface
T_surf = eos.eval_freezing_point(p0_ice, S, "Ih");
theta_surf = eos.t_to_theta(p0_ice, T_surf, S, pref);
fprintf("Surface theta: %.6f\n", theta_surf);

function poc = calc_ocean_pres(r, rho_oc, G, a, Hice, M, rho_Ih, p0_ice)
% Calculates the interior ocean pressure at radius r assuming uniform
% density rho_oc.
    r_inner = a - Hice;
    poc = rho_oc*G*(r_inner-r)/(r_inner*r) * (M - 4*pi/3*rho_Ih*(3*a*r_inner + Hice^2) ...
            - 2*pi/3*rho_oc*r_inner*(r_inner - r)*(2*r_inner + r));
    poc = poc + p0_ice;
end

% iteratively estimate rho_oc and H_oc
H_iters = 3;
rho_iters = 3;
rho_oci = 1100.0;
Hoc_i   = 60.0e+03;
trap_res = 31;

for hi = 1:H_iters
    r_coords = linspace(a - Hice - Hoc_i, a - Hice, trap_res);  % Use Hoc_i, not Hoc
    rho_p = @(p) eos.rho_pThetaS(p, theta_surf, S, pref);
    
    for rhoi = 1:rho_iters
        % Pass parameters explicitly
        pres_r = @(r) calc_ocean_pres(r, rho_oci, G, a, Hice, M, rho_Ih, p0_ice);
        integrand = zeros(size(r_coords));
        for idx = 1:length(r_coords)
            integrand(idx) = rho_p(pres_r(r_coords(idx)));
        end
        % Use Hoc_i consistently
        rho_oci = trapz(r_coords, integrand) / Hoc_i;
    end
    ac = a - Hice - Hoc_i;
    rho_cr = (1/ac^3) * (3*M/4/pi - rho_oci*((ac + Hoc_i)^3 - ac^3) - rho_Ih*(a^3 - (ac + Hoc_i)^3));
    % TODO
    delta_theta = @(p) eos.t_to_theta(p, eos.eval_freezing_point(p, S, bottom_ice_phase), S, pref) - theta_surf;
    p_guess = calc_ocean_pres(a - Hice - Hoc_i, rho_oci, G, a, Hice, M, rho_Ih, p0_ice);
    fprintf('Pressure guess %.6f yields freezing point:\n', p_guess);
    fprintf('%.6f\n', theta_surf + delta_theta(p_guess));
    
    % Use single initial guess (fzero will use secant/quasi-Newton methods)
    p_freeze = fzero(delta_theta, p_guess);
    % fprintf('Freezing pressure: %.6f\n', p_freeze);
    
    % invert pres_func to find the freezing radius
    delta_pres = @(r) calc_ocean_pres(r, rho_oci, G, a, Hice, M, rho_Ih, p0_ice) - p_freeze;
    r_guess1 = a - Hice;
    r_guess2 = 0.1 * a;
    Hoc_i = a - fzero(delta_pres, [r_guess1, r_guess2]) - Hice;
    % fprintf('Ocean depth after %d iterations: %.6f\n', mi, Hoc_i);
end
Hoc_0 = Hoc_i; rho_oc = rho_oci;
fprintf("Ocean depth: %.6f\n", Hoc_0);
fprintf("Mean ocean density: %.6f\n", rho_oc)
func_replace_string('data','rhoConst=',sprintf('rhoConst=%f,',rho_oci));

%% set up horizontal domain and legendre polynomials
nxy = nx*ny;
d2r=pi/180;
yc=-yyM+dy/2:dy:yyM; xc=[dx/2:dx:nx*dx]'-nx*dx/2;
yi=-yyM:dy:yyM; xi=[0:dx:nx*dx]'-nx*dx/2;
clat=cosd(yc);slat=sind(yc);

% Legendre polynomials in latitude
P1=slat;
P2=(3/2).*slat.^2-(1/2);
P3=(5/2).*slat.^3-(3/2).*slat;
P4=(35*slat.^4-30.*slat.^2+3)/8.0;
P5=(63.*slat.^5-70.*slat.^3+15.*slat)./8;
P6=(231.*slat.^6-315.*slat.^4+105.*slat.^2-5)/16.0;
Ps=[P1;P2;P3;P4;P5;P6];
wgt1=ny*(sind(yi(2:end))-sind(yi(1:end-1)))./(sind(yi(end))-sind(yi(1)));
wgt=repmat(wgt1,[nx,1]);

%% -- shelf ice topography & thickness
% create thickness array
ySize = size(yc);
SHice=Hice.*(yc./yc);
SHice = repmat(reshape(SHice,[1,ny]),[nx,1]);
Hunder=-SHice;
Hunder0=-(rho_Ih/rho_oc).*Hice;
fname='icetopo.bin';

fid=fopen(fname,'w','b');
fwrite(fid,Hunder,prec);
fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data.shelfice','SHELFICEtopoFile',['SHELFICEtopoFile=''',fname,''','])
fprintf(['SHELFICEtopoFile=''',fname,''','])

% ice shelf mass -- used for calculating pressure & freezing temp. 
Miceshelf=SHice.*rho_Ih;
fname='iceShelf_Mass.bin';
fid=fopen(fname,'w','b');
fwrite(fid,Miceshelf,prec);
fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data.shelfice','SHELFICEmassFile',['SHELFICEmassFile=''',fname,''','])
fprintf(['SHELFICEmassFile=''',fname,''','])

%% bathymetry & vertical profile
z_bathy = -(Hice + Hoc_0);  % total depth of model domain
if simpleBathy == 1         % bathymetry follows a sin^2 function of total amplitude dHbot
    bH2 = -(2/3)*dHbot;
    bHs = [0.,bH2,0.,0.,0.,0.,];
end
for m = 1:length(bHs)
    z_bathy = z_bathy + bHs(m)*Ps(m,:);
end
% z_bathy=repmat(reshape(z_bathy,[1,ny]),[nx,1]);
z_min = min(z_bathy);
Hliq = -(z_min + Hice);     % total "liquid depth"

ratio = 1.0 + drchfrac;
transition = [];
dz_current = dz_min;
while dz_current < dz_max
    transition = [transition; dz_current];
    dz_current = dz_current * ratio;
end
    
boundary_depth = 2.0 * sum(transition);
if boundary_depth > Hliq
    error(sprintf( ...
        ['The two boundary transition layers together require ', ...
         '%.2f m, which exceeds H = %.2f m. ', ...
         'Decrease dz_min, increase dz_max, or increase H.'], ...
    boundary_depth, Hliq));
end
    
interior_depth = Hliq - boundary_depth;
n_interior = round(interior_depth / dz_max);

% Adjust dz_max to fill interior depth exactly
if n_interior > 0
    dz_interior = interior_depth / n_interior;
else
    % No interior cells — rescale transition layer uniformly to fill H
    scale = Hliq / boundary_depth;
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
dh = transpose([Hice; dz]);
fprintf("dh profile (km): %.6f\n", dh/1.0e+3);
hf = [0,cumsum(dh)];
hc=(hf(2:end)+hf(1:end-1))./2;  % centered depth coordinates
rc=a-hc;
nr=length(hc);
fprintf("Domain size in z %i\n:",nr)

zbathy_twod = repmat(z_bathy, nx, 1);
zbathy_twod(:,1)=zeros(nx,1);
zbathy_twod(:,ny)=zeros(nx,1);
file_name=['my_bathy_',sprintf('%d',-z_min/1.0e+3),'km.bin'] ;
fid=fopen(file_name,'w','b');
fwrite(fid,zbathy_twod,prec); fclose(fid);
fprintf(['\nwrite file: ', file_name])
func_replace_string('data','bathyFile',['bathyFile=''',file_name,''','])
fprintf(['bathyFile=''',file_name,''','])


%% gravity profile

g0 = G*M/a^2;
g_top = p0_ice/(rho_Ih*dh(1));
g = [g_top, 4*G*pi/3.*((rho_cr - rho_oc).*ac^3./rc(2:nr).^2 + rho_oc.*rc(2:nr))];
grav_fac = g/g0;
% fprintf("gravity profile (m/s^2): %.6f\n", g);

func_replace_string('data','gravity=',sprintf('gravity=%f,',g0));
fname=strcat('gravity_', satellite_name,'.bin');
fid=fopen(fname,'w','b');
fwrite(fid,g,prec);
fclose(fid);
fprintf(['\nwrite file: ',fname]);
func_replace_string('data','gravityFile',['gravityFile=''',fname,''','])
fprintf(['gravityFile=''',fname,''','])

%% initial temperature profile
if Tini_profile == 0
    Tprof = theta_surf*ones(nr);
elseif Tini_profile == 1
    Tprof = theta_surf*ones([1,nr]);
    hbathy_max = - max(z_bathy);
    mask=hf(2:end) > hbathy_max;
    bathy_rs = rc(mask);
    bathy_ps = arrayfun(@(r) calc_ocean_pres(r, rho_oci, G, a, Hice, M, rho_Ih, p0_ice), bathy_rs);
    bathy_thetas = arrayfun(@(p) eos.t_to_theta(p, eos.eval_freezing_point(p, S, bottom_ice_phase), S, pref), bathy_ps);
    Tprof(mask) = bathy_thetas;
    Tprof(~mask) = min(bathy_thetas);
elseif Tini_profile == 2
    Tini = restart_data.theta;
    Tini = permute(Tini, ndims(Tini):-1:1);
    Tini = flip(Tini,3);
    sfx='pickup.bin';
    fTname=['theta_ini.',sfx]
end

if Tini_profile < 2
    % fprintf("Vertical temperature profile (K): %.6f\n", Tprof);
    Tini=repmat(Tprof, nxy, 1);
    % add random noise
    delTemp = 1.0e-3;
    dTn=rand(size(Tini));
    Tini=Tini + delTemp*dTn - 273.15;
    sfx=sprintf('%i%s',delTemp*1.0e3,'mK.bin');
    fTname=['theta_ini.',sfx];
end 
fid=fopen(fTname,'w','b'); fwrite(fid,Tini,prec); fclose(fid);
fprintf(['\nwrite file: ',fTname]);
func_replace_string('data','hydrogThetaFile',['hydrogThetaFile=''',fTname,''','])
fprintf(['hydrogThetaFile=''',fTname,''','])

%% initial salt field
if Sini_profile == 2
    Sini = restart_data.S;
    Sini = permute(Sini, ndims(Sini):-1:1);
    Sini = flip(Sini,3);
    sfx='pickup.bin';
    fSname=['salt_ini.',sfx]
    fid=fopen(fSname,'w','b'); fwrite(fid,Sini,prec); fclose(fid);
    fprintf(['\nwrite file: ',fSname]);
    func_replace_string('data','hydrogSaltFile',['hydrogSaltFile=''',fSname,''','])
    fprintf(['hydrogSaltFile=''',fSname,''','])
end

%% initial velocity fields
if init_Vel == 1
    Uini = restart_data.U; Vini = restart_data.V;
    Uini = flip(permute(Uini, ndims(Uini):-1:1),3);
    Vini = flip(permute(Vini, ndims(Vini):-1:1),3);
    sfx='pickup.bin';
    fUname=['uVel_ini.',sfx]; fVname=['vVel_ini.',sfx];
    fid=fopen(fUname,'w','b');fwrite(fid,Uini,prec);fclose(fid);
    fid=fopen(fVname,'w','b');fwrite(fid,Vini,prec);fclose(fid);
    func_replace_string('data','uVelInitFile',['uVelInitFile=''',fUname,''','])
    func_replace_string('data','vVelInitFile',['vVelInitFile=''',fVname,''','])
    fprintf(['uVelInitFile=''',fUname,''','])
    fprintf(['vVelInitFile=''',fVname,''','])
end

%% bottom heating profile
if prescribe_Qbot > 0
    func_replace_string('data.seafrz','sfz_QbotType',['sfz_QbotType=1,'])
    if prescribe_Qbot == 1
        Qbot = 1.0e-2;
        Qbot = repmat(Qbot,[nx,ny]);
        fname = 'Qbot_uniform.bin';
    end
    fid = fopen(fname,'w','b');
    fwrite(fid,Qbot,prec);
    fclose(fid);
    func_replace_string('data.seafrz','sfz_QbotFile',['sfz_QbotFile=''',fname,''','])
    fprintf(['sfz_QbotFile=''',fname,''','])
end

end

