% equation of state calculations for gendata.m

classdef eosMethods
    properties
        Gcoeffs
        normP
        normS
        normT
        Ih
        III
        V
        VI
        adtg_coeffs
    end
    methods
        function obj = eosMethods()
            gibbs_data = load('eosdata/seawater_gibbs_coeffs.mat')
            obj.Gcoeffs = gibbs_data.Gcoeffs;
            obj.normP   = gibbs_data.normP;
            obj.normS   = gibbs_data.normS;
            obj.normT   = gibbs_data.normT;
            adtg_data = load('eosdata/adtg_coeffs.mat');
            obj.adtg_coeffs = adtg_data.adtg_coeffs; 
            fp_data = load('eosdata/melting_curve_linear_params.mat');
            obj.Ih = fp_data.Ih;
            obj.III = fp_data.III;
            obj.V = fp_data.V;
            obj.VI = fp_data.VI;
        end

        %% eos equations
        
        function rho = rho_pTS(obj, p, T, S)
        % Evaluates the density in kg/m^3 at pressure p [Pa], temperature T [K],
        % salinity S [g/kg].
            vspec = obj.eval_tricubic(p, T, S, obj.normP(1), obj.normP(2), obj.normT(1), obj.normT(2), obj.normS(1), obj.normS(2), obj.Gcoeffs, 1, 0, 0);
            rho = 1/vspec;
        end
        
        function fp = eval_freezing_point(obj, p, S, ice_type)
            if ice_type == "Ih"
                fp = obj.Ih(1)*p + obj.Ih(2)*S + obj.Ih(3);
            elseif ice_type == "III"
                fp = obj.III(1)*p + obj.III(2)*S + obj.III(3);
            elseif ice_type == "V"
                fp = obj.V(1)*p + obj.V(2)*S + obj.V(3); 
            elseif ice_type == "VI"
                fp = obj.VI(1)*p + obj.VI(2)*S + obj.VI(3);
            end
        end
        
        %% potential temperature
        function theta = t_to_theta(obj, p, T, S, pref)
            adtg_coeffs = obj.adtg_coeffs;
            gamma = adtg_coeffs(1)+S*adtg_coeffs(2) ...
                +0.5*(p+pref)*(adtg_coeffs(3)+S*adtg_coeffs(4)) ...
                +T*(adtg_coeffs(5)+S*adtg_coeffs(6)+T*adtg_coeffs(7)+p*adtg_coeffs(8));
            theta = T - (p - pref)*gamma;
        end
        
        function t = theta_to_t(obj, p, theta, S, pref)
            t = obj.t_to_theta(pref, theta, S, p);
        end

        function rho = rho_pThetaS(obj, p, theta, S, pref)
            T = obj.theta_to_t(p, theta, S, pref);
            rho = obj.rho_pTS(p, T, S);
        end
    end
    methods(Static)
        % tricubic polynomial evaluator 
        function result = eval_tricubic(xq, yq, zq, xc, xh, yc, yh, zc, zh, coeffs, mx, my, mz)
        % EVAL_TRICUBIC  Evaluate a fitted tricubic polynomial and its partial derivatives.
        %
        % USAGE:
        %   result = eval_tricubic(xq, yq, zq, xc, xh, yc, yh, zc, zh, coeffs)
        %   result = eval_tricubic(xq, yq, zq, xc, xh, yc, yh, zc, zh, coeffs, mx, my, mz)
        %
        % INPUTS:
        %   xq, yq, zq  -- query points (scalars or broadcast-compatible arrays)
        %   xc, xh      -- x normalisation: centre and half-range
        %   yc, yh      -- y normalisation: centre and half-range
        %   zc, zh      -- z normalisation: centre and half-range
        %   coeffs      -- (4,4,4) array of fitted coefficients c(p,q,r),
        %                  where p,q,r are 1-indexed in Matlab (i.e. p=1..4)
        %   mx, my, mz  -- (optional) derivative orders w.r.t. x, y, z; default 0.
        %                  Each must be in {0,1,2,3}.
        %
        % OUTPUT:
        %   result      -- d^(mx+my+mz) f / dx^mx dy^my dz^mz evaluated at (xq,yq,zq).
        %                  Same shape as the broadcast of xq, yq, zq.
        %
        % NOTES:
        %   Polynomial model:
        %     f(x,y,z) = sum_{p,q,r=0}^{3} coeffs(p+1,q+1,r+1) * xn^p * yn^q * zn^r
        %   where xn = (x-xc)/xh, yn = (y-yc)/yh, zn = (z-zc)/zh.
        %
        %   Differentiation replaces xn^p with its m-th derivative:
        %     d^m/dx^m xn^p = (1/xh)^m * [p!/(p-m)!] * xn^(p-m)   for p >= m
        %                   = 0                                       for p <  m
        
            % --- Default derivative orders ---
            if nargin < 11, mx = 0; end
            if nargin < 12, my = 0; end
            if nargin < 13, mz = 0; end
        
            if any([mx my mz] < 0) || any([mx my mz] > 3)
                error('Derivative orders mx, my, mz must each be in {0,1,2,3}.');
            end
        
            % --- Falling-factorial table: ffac(p,m) = p!/(p-m)!, p=1..4, m=1..4 ---
            % Row index is p (1-indexed), column index is m (1-indexed).
            % Entry is 0 when m > p (i.e. when the derivative annihilates the monomial).
            ffac = [1  0  0  0;   % p=0
                    1  1  0  0;   % p=1
                    1  2  2  0;   % p=2
                    1  3  6  6];  % p=3
            % ffac(p+1, m+1) gives the falling factorial for 0-indexed p,m.
        
            % --- Normalise coordinates ---
            xn = (xq - xc) / xh;
            yn = (yq - yc) / yh;
            zn = (zq - zc) / zh;
        
            % --- Chain-rule scale factors ---
            scale = (1/xh)^mx * (1/yh)^my * (1/zh)^mz;
        
            % --- Derivative power vectors ---
            % dpx(p+1) = [p!/(p-m)!] * xn^(p-m)  for p >= m, else 0.
            % Shape matches xn so this works for scalar or array inputs.
            dpx = eosMethods.deriv_powers(xn, mx, ffac);   % (4, ...) where ... is size(xn)
            dpy = eosMethods.deriv_powers(yn, my, ffac);
            dpz = eosMethods.deriv_powers(zn, mz, ffac);
        
            % --- Evaluate sum_{p,q,r} coeffs(p,q,r) * dpx(p) * dpy(q) * dpz(r) ---
            result = zeros(size(xn));
            for r = 1:4
                for q = 1:4
                    for p = 1:4
                        result = result + coeffs(p,q,r) .* dpx(p,:) .* dpy(q,:) .* dpz(r,:);
                    end
                end
            end
        
            result = scale .* result;
        
        end
        
        
        % -------------------------------------------------------------------------
        function dp = deriv_powers(un, m, ffac)
        % DERIV_POWERS  Compute differentiated monomial values for a 1-D coordinate.
        %
        % Returns dp of size (4, numel(un)) where:
        %   dp(p+1, :) = ffac(p,m) * un.^(p-m)   for p >= m
        %   dp(p+1, :) = 0                         for p <  m
        %
        % un   -- normalised coordinate (any shape; treated as a flat vector)
        % m    -- derivative order (0..3)
        % ffac -- (4,4) falling-factorial table (0-indexed p,m stored at p+1,m+1)
        
            un_flat = un(:).';          % (1, N) row vector for clean broadcasting
            dp = zeros(4, numel(un));
        
            for p = 0:3
                if p >= m
                    dp(p+1, :) = ffac(p+1, m+1) .* un_flat .^ (p - m);
                end
                % p < m case is already zero from initialisation
            end
        
        end
        % -------------------------------------------------------------------------
    end
end
