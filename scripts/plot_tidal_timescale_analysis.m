function plot_tidal_timescale_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% save_flag: bool

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));

% Stability criteria
a_out_stability_Vynatheya_circular   = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
% # 0.1 Z_SMC
% idx_0_1_Z_SMC = [47, 1784, 2132, 2504]
% age_yr_0_1_Z_SMC = [5.329669e+03, 5.120677e+06, 5.307573e+06, 5316369.558564772]
% apsidal_constant_k2_0_1_Z_SMC = [0.023533, 0.008881, 0.004897, 0.0007136587012452427]
% period_days_0_1_Z_SMC = [1.099654, 1.450562, 1.668519, 1.6824947859764587]
% radius_Rsun_0_1_Z_SMC = [7.577242, 2.401224, 1.939817, 0.8534730543749688]
% mass_conv_core_Msun_0_1_Z_SMC = [37.805430, 41.412994, 38.487511, 0.0]
% mass_Msun_0_1_Z_SMC = [54.999836, 43.976249, 43.793192, 43.79319230382309]
% Choose
filename            = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR_Tides.mat';
title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.00035/ \rm{at\ ZAMS}$'];
plot_label_temp_png = '../plots/png/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_';
plot_label_temp_pdf = '../plots/pdf/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_';    
mass_Msun               = 55;
radius_Rsun             = round(7.577242,2);
orbital_period_days     = 1.1; 
% Double check the following values with MESA
mass_conv_core_Msun     = round(37.805430,2);
radius_conv_core_Rsun   = 3.5; % Approximate from MESA model, without going into profile
luminosity_Lsun         = 10^5.632404; % Pending value

% Load
if debug_flag==true
    M                               = load(filename)
else
    M                               = load(filename);
end

e_max       = M.eccs;
m3          = M.m3;
Pout        = M.p2;
tau_sec     = M.tau_sec.*AstroConstants.s_to_yr;

% Calculate extra values
orbital_period_year = orbital_period_days./AstroConstants.yr_to_d;
separation_inner_AU = separation_in_AU(orbital_period_year,mass_Msun+mass_Msun);

tau_circ_dyn_yr = calculate_timescale_dynamical_tides(mass_Msun, radius_Rsun, separation_inner_AU, mass_Msun)
tau_circ_equ_yr = calculate_timescale_equilibrium_tides(mass_Msun, mass_conv_core_Msun, radius_Rsun, radius_conv_core_Rsun, luminosity_Lsun, separation_inner_AU, mass_Msun)

% Create 2D meshgrid
[X, Y]                          = meshgrid(m3, Pout);

% Calculate eta=1
mock_mass = linspace(1,100,1001);

q_out                       = mock_mass./(mass_Msun+mass_Msun);
crit_stability_Vynatheya_a_AU      = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr  = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d   = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;
crit_stability_P_orb_d = crit_stability_Vynatheya_P_orb_d;

% Filtering
idx_dyn_tides = find(tau_sec>=tau_circ_dyn_yr);
idx_equ_tides = find(tau_sec>=tau_circ_equ_yr);

tau_sec(idx_dyn_tides) = nan;
tau_sec(idx_equ_tides) = nan;

% e_lim_val = 1-(2*radius_Rsun/(separation_inner_AU.*AstroConstants.AU_to_Rsun));
% e_lim_tol = 0.01;
% idx_of_eccentricity = find((e_max'>=e_lim_val-e_lim_tol)&(e_max'<=e_lim_val+e_lim_tol)&(Pout'>=5));
% 
% min_eccentricity = 0.001;
% index_of_non_induced_eccentricity = find(e_max<min_eccentricity);
% cos_inc(index_of_non_induced_eccentricity) = nan;
% e_max(index_of_non_induced_eccentricity) = nan;
% 
% min_val_colorbar    = log10(min(min(min(min(eps_SA)),min(min(eps_GR))),min(min(eps_Tide))));
% max_val_colorbar    = log10(max(max(max(max(eps_SA)),max(max(eps_GR))),max(max(eps_Tide))));
% extreme_val_colorbar= max(abs(min_val_colorbar),abs(max_val_colorbar));

% Print values
if debug_flag==true
    title_string
    fprintf('Mass_1 = Mass 2 = %f Msun',mass_Msun)
    fprintf('\n')
    fprintf('a_{inner} = %f au', separation_inner_AU)
    fprintf('\n')
    fprintf('P_{orb} = %f d', orbital_period_days)
    fprintf('\n')    
end

% PLOT
list_plot = {'\tau_{sec}'};
[indx_plot, tf_plot] = listdlg('ListString',list_plot);

if indx_plot==1
    plot_label_png = strcat(plot_label_temp_png,'tau_sec.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'tau_sec.pdf');
else
    warning("Odd choice.")
end  


sz = 145.0;
sz2 = 20;
fs=18;

lines6      = [    0.3010    0.7450    0.9330];
dark_grey   = 0.3.*[1 1 1];
light_grey  = 0.6.*[1 1 1];    
grey   = 0.5.*[1 1 1];

xLims       = [1 100];
xLimsLabels = {'1','10','100'};
yLims       = [1 1000];
yLimsLabels = {'1','10','100','1000'};

extreme_val_tau = max(abs(log10(min(min(tau_sec)))),abs(log10(max(max(tau_sec)))));

leveler = 10;

clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')

hold on        
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink(1000)))
cbar.Label.Interpreter = 'latex';
cbar.XTick = [-4:2:4];

if indx_plot == 1
    % secular timescale
    cbar.Label.String = '$\log_{10} (\tau_{\rm{sec}}/\rm{yr})$';
    mm_tau_sec=mesh(X, Y, log10(tau_sec'),'HandleVisibility','off');
    mm_tau_sec.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_tau,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_tau,'HandleVisibility','off')
    empty_region = fill3([1 100 100 1],[10 10 1000 1000 ], -2.*[1 1 1 1],'k');    
    colormap(flip(slanCM('vik')))  
else    
    warning("Odd choice.")    
end

if annotations_flag
    empty_region.EdgeColor = 'none';
    empty_region.FaceColor = light_grey;
    empty_region.HandleVisibility = 'off';
    text(1.2,600,leveler,'$\tau_{\rm{sec}} > \tau_{\rm{tides}}$','Color','w','Fontsize',fs)
end

dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],dark_grey);
unstableRegion.EdgeColor = 'none';
unstableRegion.FaceColor = dark_grey;
unstableRegion.HandleVisibility = 'off';
text(1.2,1.75,leveler,'Dynamically Unstable','Color','w','Fontsize',fs)

ylabel('$P_{\rm{out}}/\rm{d}$')
xlabel('$m_3/M_{\odot}$')
pbaspect([2 1 1])
ax1 = gca;
ax1.Layer = 'top';
ax1.FontSize = fs;
ax1.XLim = xLims;
ax1.YLim = yLims;
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.Layer = 'top';
ax1.XTickLabel = xLimsLabels;
ax1.YTickLabel = yLimsLabels;
ax1.TickLabelInterpreter = 'latex';
set(cbar,'TickLabelInterpreter','latex')
box on
  
% Save
if save_flag
    print(gcf,plot_label_png,'-dpng','-r250');
    print(gcf,plot_label_pdf,'-dpdf');    
end

end

function tau_circ_dyn_yr = calculate_timescale_dynamical_tides(M_Msun, R_Rsun, a_AU, M_comp)
% Based on Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
q2 = M_comp./M_Msun;
R_AU = R_Rsun.*AstroConstants.Rsun_to_AU;

E_2 = calculate_E_2(M_Msun);

inv_tau_dyn_circ = (21.0./2).*sqrt((AstroConstants.G_AUMsunyr.*M_Msun./(R_AU.^3))).*q2.*((1+q2).^(11.0/6)).*E_2.*((R_AU/a_AU).^(21.0/2)); % Eq. 41 of Hurley, Pols & Tout (2002)
tau_circ_dyn_yr = 1./inv_tau_dyn_circ;
end

function E_2 = calculate_E_2(M_Msun)
% Based on Zahn (1975)
% As implemented in Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
E_2 = 1.592*(10^-9).*(M_Msun.^2.84); % Eq. 43 of Hurley, Pols & Tout (2002)
end

function tau_circ_eq_yr = calculate_timescale_equilibrium_tides(M_Msun, Mconv_Msun, R_Rsun, Rconv_Rsun, L_Lsun, a_AU, M_comp)
% Based on Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
q2 = M_comp./M_Msun;
R_AU = R_Rsun.*AstroConstants.Rsun_to_AU;

t_conv_yr = 0.4311.*(Mconv_Msun.*Rconv_Rsun.*(R_Rsun-0.5.*Rconv_Rsun)./(3.0.*L_Lsun)).^(1.0/3);

% Following Preece et al. (2022)
% https://ui.adsabs.harvard.edu/abs/2022ApJ...933...25P/abstract
% k_over_T = (2.0./21.0).*(f_conv./t_conv_yr).*(Mconv_Msun/M_Msun);
chosen_metallicity = 0.001; 
a = (-0.12*log10(chosen_metallicity))+6.91;
b = (0.23*log10(chosen_metallicity))-0.5;
c = (-0.28*log10(chosen_metallicity))+0.07;
k_over_T_c = ((Rconv_Rsun./R_Rsun).^a).*((Mconv_Msun./M_Msun).^b).*c./t_conv_yr;
inv_tau_eq_circ_yr = (21.0./2).*(k_over_T_c).*q2.*(1+q2).*((R_AU/a_AU).^8.0);
tau_circ_eq_yr = 1./inv_tau_eq_circ_yr;

end