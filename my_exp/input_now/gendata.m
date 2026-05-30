% Updated version of Wanying Kang's gendata.m code for planetary structure and MITgcm initialization.
function gendata()
%% load packages
MITROOT='/home/cgliddon/'; % MITROOT should point to the directory where MITgcm is installed XX.
addpath([MITROOT,'/MITgcm/utils/matlab/'])
eos = eosMethods(); % defined in the file eosMethods.m, in same directory as gendata.m
%% Constants and free parameters
G=6.674e-11;           % Newton's gravitational constant, N m^2 kg^-2
rho_Ih=917.;           % ice-Ih density, kg m^-3
k_Ih=651;          % ice-Ih base thermal conductivity, W/m
kappa_wh=1.0e-3;       % water thermal diffusivity
kappa_wr=1.0e-3;
nu_wh=1.;              % water (eddy) viscosity, m^2 s^-1
nu_wr=1.;
alpha_T=4.5e-4;        % water thermal expansion coefficient
cpw=4.184e+03;         % water heat capacity
M=1.4819e+23;          % planetary mass, kg
a=2634.1e+03;          % planetary radius, m
Qflx_avg=1.0e-02;      % mean outward heat flux, W m^-2
Tsurf=110;             % surface temperature, K
S=10.0;                % salinity, g kg^-1
pref=1.0e+05;          % reference pressure for PT-T interconversion
bottom_ice_phase="V"; % reference bottom ice phase
delta_ht=1.0e3;        % topographic variation scale for top boundary
delta_hb=1.0e3;        % ditto, but at the bottom
boundaryTopo=0;        % 0: "pseudo" boundary topography from prescribed temperature profiles; 1: "real" topography
simpleBathy=1;         % 0: all Legendre coefficients need to be specified manually; 1: bathymetry is just -P2(sin(phi))
Tini_profile=1;        % 0: isentropic; 1: superadiabatic following Gastine et al. scaling; ; 2: read from file
Sini_profile=0;		   % 0: isohaline; 2: read from file
init_Vel=0;		       % 0: does nothing; 1: read from file
prescribe_Qbot=1;	   % 0: does nothing; 1: generates and writes out a bottom heating pattern
satellite_name='myMoon';    
restart_datfile='';     
if not(strcmp(restart_datfile, ''))
    restart_data=load(restart_datfile);
end

nx=96;ny=96;        % grid size in x and y
dx=1.75;dy=1.75; yyM=84.; % grid resolution in x and y; maximum latitude (+/-)
dz_max = 2.5e+03; dz_min = 8.0e+02;   % max and min grid spacing for vertical levels
drchfrac = 0.10;                      % how much do we allow the vertical spacing to change by for adjacent cells?

prec='real*4';                        % data save-out precision

%% start replacing parameters in data file
eosType='SEAFRZ';
func_replace_string('data','rSphere=',sprintf('rSphere=%f,',a))
func_replace_string('data','eosType',sprintf('eosType=%s,',eosType))
%func_replace_string('data','viscAh',sprintf('viscAh=%f,',nu_wh))
%func_replace_string('data','viscAr',sprintf('viscAr=%f,',nu_wr))
%func_replace_string('data','diffKhT',sprintf('diffKhT=%f,',kappa_wh))
%func_replace_string('data','diffKrT',sprintf('diffKrT=%f,',kappa_wr))
if boundaryTopo==0
    simpleBathy=1;      % override
    func_replace_string('data.hpimc','hpimc_prescribeTb',sprintf('hpimc_prescribeTb=.TRUE.,'))
else
    func_replace_string('data.hpimc','hpimc_prescribeTb',sprintf('hpimc_prescribeTb=.FALSE.,'))
end


%% Part 1: calculate balanced 1D internal structure
power=Qflx_avg*4*pi*a^2;    % total power output by radiogenic heating

Hi0=50.0e+3; Tt0=270.0;
function F = myEquations(res)
    Hi=res(1);Ti=res(2);
% we multiply by 1e6 because otherwise the residual is very small by construction, which causes problems with fsolve's tolerance    
    F(1)=(k_Ih/Hi*log(Ti/Tsurf) - Qflx_avg)*1.0e6;   
    p1=(rho_Ih*G*Hi)/(a^2-a*Hi)*(M+2*pi*rho_Ih*a^2*Hi-2*pi*rho_Ih*a*Hi^2/3);
%    fprintf('p1=%.6f\n',p1);
    F(2)=Ti-eos.eval_freezing_point(p1,S,"Ih");
end 
x0 = [Hi0,Tt0];
sol = fsolve(@myEquations, x0);
Hice0=sol(1);Tt=sol(2);
p1=(rho_Ih*G*Hice0)/(a^2-a*Hice0)*(M+2*pi*rho_Ih*a^2*Hice0-2*pi*rho_Ih*a*Hice0^2/3);
theta_t=eos.t_to_theta(p1,Tt,S,pref);
fprintf("Ice shell thickness: %.6f\n", Hice0)
fprintf("Surface theta: %.6f\n", theta_t)

% gravity at the top of the ocean
g_top = p1/(rho_Ih*Hice0);

% define function for nonlinear solving
function poc = calc_ocean_pres(r, rho_oc, G, a, Hice0, M, rho_Ih, p1)
% Calculates the interior ocean pressure at radius r assuming uniform
% density rho_oc.
    r_inner = a - Hice0;
    poc = rho_oc*G*(r_inner-r)/(r_inner*r) * (M - 4*pi/3*rho_Ih*(3*a*r_inner + Hice0^2) ...
            - 2*pi/3*rho_oc*r_inner*(r_inner - r)*(2*r_inner + r));
    poc = poc + p1;
end

if Tini_profile==0 | Tini_profile==1
    % iteratively estimate rho_oc and H_oc
    H_iters = 5;
    rho_iters = 3;
    rho_oci = 1100.0;
    Hoc_i   = 60.0e+03;
    deltav_theta = 0.02;
    trap_res = 31;

    for hi = 1:H_iters
        r_coords = linspace(a - Hice0 - Hoc_i, a - Hice0, trap_res);  % Use Hoc_i, not Hoc
        rho_p = @(p) eos.rho_pThetaS(p, theta_t, S, pref);
        
        for rhoi = 1:rho_iters
            % Pass parameters explicitly
            pres_r = @(r) calc_ocean_pres(r, rho_oci, G, a, Hice0, M, rho_Ih, p1);
            integrand = zeros(size(r_coords));
            for idx = 1:length(r_coords)
                integrand(idx) = rho_p(pres_r(r_coords(idx)));
            end
            % Use Hoc_i consistently
            rho_oci = trapz(r_coords, integrand) / Hoc_i;
        end
        ac = a - Hice0 - Hoc_i;
        rho_cr = (1/ac^3) * (3*M/4/pi - rho_oci*((ac + Hoc_i)^3 - ac^3) - rho_Ih*(a^3 - (ac + Hoc_i)^3));
        % TODO
        delta_theta = @(p) eos.t_to_theta(p, eos.eval_freezing_point(p, S, bottom_ice_phase), S, pref) - deltav_theta - theta_t;
        p_guess = calc_ocean_pres(a - Hice0 - Hoc_i, rho_oci, G, a, Hice0, M, rho_Ih, p1);
        fprintf('Pressure guess %.6f yields freezing point:\n', p_guess);
        fprintf('%.6f\n', theta_t + deltav_theta + delta_theta(p_guess));
        
        % Use single initial guess (fzero will use secant/quasi-Newton methods)
        p_freeze = fzero(delta_theta, p_guess);
        % fprintf('Freezing pressure: %.6f\n', p_freeze);
        
        % invert pres_func to find the freezing radius
        delta_pres = @(r) calc_ocean_pres(r, rho_oci, G, a, Hice0, M, rho_Ih, p1) - p_freeze;
        r_guess1 = a - Hice0;
        r_guess2 = 0.1 * a;
        Hoc_i = a - fzero(delta_pres, [r_guess1, r_guess2]) - Hice0;
        % fprintf('Ocean depth after %d iterations: %.6f\n', mi, Hoc_i);
        if Tini_profile == 1
	        deltav_theta=sqrt(nu_wr*(ac/(a-Hice0))*power/(alpha_T*rho_oci*cpw*g_top*Hoc_i^2));
        end
        fprintf('Superadiabatic temperature gradient (K): %.6f\n',deltav_theta)
    end
% else if Tini_profile == 1
end
Hoc0 = Hoc_i; rho0 = rho_oci; p2=p_freeze;
fprintf("Ocean depth: %.6f\n", Hoc0);
fprintf("Mean ocean density: %.6f\n", rho0)
func_replace_string('data','rhoConst=',sprintf('rhoConst=%f,',rho0));

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

%% shelf ice topography & thickness
if boundaryTopo==0
    SHice=Hice0.*(yc./yc);
    SHice = repmat(reshape(SHice,[1,ny]),[nx,1]);
    Hunder=-SHice;
end
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
if boundaryTopo==0
    % no topography, just 
    zmin=-(Hice0+Hoc0);
    Hliq=Hoc0;
    zbathy=zmin;
end

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
dh = transpose([Hice0; dz]);
fprintf("dh profile (km): %.6f\n", dh/1.0e+3);
%func_replace_string('data','delR',['delR=',sprintf('%f,',dh),',']);
hf = [0,cumsum(dh)];
hc=(hf(2:end)+hf(1:end-1))./2;  % centered depth coordinates
rc=a-hc;
nr=length(hc);
fprintf("Domain size in z %i:\n",nr)

zbathy_twod = repmat(zbathy,[nx,ny]);
zbathy_twod(:,1)=zeros(nx,1);
zbathy_twod(:,ny)=zeros(nx,1);
file_name=['my_bathy_',sprintf('%d',-zmin/1.0e+3),'km.bin'] ;
fid=fopen(file_name,'w','b');
fwrite(fid,zbathy_twod,prec); fclose(fid);
fprintf(['\nwrite file: ', file_name])
func_replace_string('data','bathyFile',['bathyFile=''',file_name,''','])
fprintf(['bathyFile=''',file_name,''','])


%% gravity profile

g0 = G*M/a^2;
g_top = p1/(rho_Ih*dh(1));
g = [g_top, 4*G*pi/3.*((rho_cr - rho0).*ac^3./rc(2:nr).^2 + rho0.*rc(2:nr))];
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
if Tini_profile == 0 | Tini_profile == 1
    Tprof = theta_t*ones(1, nr);
%    deltav_theta=0;
    Tprof(2:end)=Tprof(2:end)+deltav_theta*(hc(2:end)-mean(hc(2:end)))/(hf(end)-hf(2));
elseif Tini_profile == 3
    Tprof = theta_t*ones([1,nr]);
    hbathy_max = - max(zbathy);
    mask=hf(2:end) > hbathy_max;
    bathy_rs = rc(mask);
    bathy_ps = arrayfun(@(r) calc_ocean_pres(r, rho_oci, G, a, Hice0, M, rho_Ih, p1), bathy_rs);
    bathy_thetas = arrayfun(@(p) eos.t_to_theta(p, eos.eval_freezing_point(p, S, bottom_ice_phase), S, pref), bathy_ps);
    Tprof(mask) = bathy_thetas;
    Tprof(~mask) = min(bathy_thetas);
    deltav_theta=Tprof(end)-Tprof(1);
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
    func_replace_string('data.hpimc','hpimc_QbotType','hpimc_QbotType=1,')
    if prescribe_Qbot == 1
        Qbot = repmat(Qflx_avg,[nx,ny]);
        fname='Qbot_uniform.bin';
        fid = fopen(fname,'w','b');
        fwrite(fid,Qbot,prec);
        fclose(fid);
        func_replace_string('data.hpimc','hpimc_QbotFile',['hpimc_QbotFile=''',fname,''','])
        fprintf(['hpimc_QbotFile=''',fname,''','])
    end
end

%% bottom temperature profile
if boundaryTopo == 0
    % top boundary
    rho1=eos.rho_pThetaS(p1, theta_t, S, pref);
    topo1=delta_ht.*P2;   % multiply by 2nd Legendre polynomial
    deltap1=rho1*g_top*topo1;
    T1=arrayfun(@(p) eos.eval_freezing_point(p, S, 'Ih'), deltap1+p1);
    theta1=arrayfun(@(T) eos.t_to_theta(p1, T, S, pref), T1);
    theta1=repmat(theta1,[nx,1])-273.15;
    % bottom boundary
    theta_b=theta_t+deltav_theta;
    rho2=eos.rho_pThetaS(p2,theta_b,S,pref);
    topo2=-delta_hb.*P2;
    deltap2=rho2*g(end)*topo2;
    T2=arrayfun(@(p) eos.eval_freezing_point(p, S, bottom_ice_phase), deltap2+p2);
    theta2=arrayfun(@(T) eos.t_to_theta(p2, T, S, pref), T2);
    theta2=repmat(theta2,[nx,1])-273.15;
    fprintf("Ocean top boundary layer PT: %.6f\n", mean(theta1));
    fprintf("Ocean bottom boundary layer PT: %.6f\n", mean(theta2));

    % save out
    fname1='thetaBL_top.bin';fname2='thetaBL_bot.bin';
    fid=fopen(fname1,'w','b'); fwrite(fid,theta1,prec); fclose(fid);
    fid=fopen(fname2,'w','b'); fwrite(fid,theta2,prec); fclose(fid);
    func_replace_string('data.hpimc', 'hpimc_TbtopFile', ['hpimc_TbtopFile=''',fname1,''',']);
    func_replace_string('data.hpimc', 'hpimc_TbbotFile', ['hpimc_TbbotFile=''',fname2,''',']);
    fprintf(['hpimc_TbtopFile=''',fname1,''',']);
    fprintf(['hpimc_TbbotFile=''',fname2,''',']);
end


end

