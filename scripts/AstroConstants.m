classdef AstroConstants
   properties (Constant)
       c_cgs        = 2.99792458*10^10;     % speed of light in cm/s
       c_mks        = 2.99792458e8          % speed of light in m/s 
       G_mks        = 6.67259e-11           % Gravitational constant in m/kgs
       G_AUMsunyr   = 4*pi*pi               % Gravitational constant AU^3 yr^-2 M^-1  
       Msun_kg      = 1.9891e30             % Msun in kg
       Msun_g       = 1.9891e+33;           % Msun g
       Rsun_m       = 6.957e+8              % Rsun in m
       AU_to_m      = 1.496e+11             % AU in m
       AU_to_Rsun   = 215.032               % AU in Rsun
       Rsun_to_AU   = 1./215.032;        % Rsun in AU
       s_to_yr      = 3.171e-8;             % seconds to year 
       yr_to_d      = 365.25;               % yr in days
       yr_to_Myr    = 1e-6;                 % yr in Myr
       km_to_cm     = 1e+5;                 % km to cm
       t_Hubble_MYR = 13.968*1e3;           % Hubble time in Myr
   end
end