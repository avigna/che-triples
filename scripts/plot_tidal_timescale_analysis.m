function plot_tidal_timescale_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gómez
% debug_flag: bool
% annotations_flag: bool
% save_flag: bool

% FUNCTIONS
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));

% Stability criteria
a_out_stability_Vynatheya_circular   = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
filename            = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR_Tides.mat';
title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.00035/ \rm{at\ ZAMS}$'];
plot_label_temp_png = '../plots/png/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_low_Z_';
plot_label_temp_pdf = '../plots/pdf/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_low_Z_';    
mass_Msun               = 55;
radius_Rsun             = round(7.577242,2);
orbital_period_days     = 1.1; 
mass_conv_core_Msun     = round(37.805430,2);
radius_conv_core_Rsun   = 3.5; % Approximate from MESA model, without going into profile
k_2 = 0.024;

% Load
if debug_flag==true
    M                               = load(filename)
else
    M                               = load(filename);
end

eps_Tide    = M.eps_Tide;
m3          = M.m3;
Pout        = M.p2;
tau_ZLK     = M.tau_sec.*AstroConstants.s_to_yr.*3.0./4; % From Evgeni's script; extra factor to correct

% Calculate extra values
orbital_period_year = orbital_period_days./AstroConstants.yr_to_d;
separation_inner_AU = separation_in_AU(orbital_period_year,mass_Msun+mass_Msun);

tau_dyn_sync_yr = calculate_timescale_dynamical_tides_sync(mass_Msun, mass_conv_core_Msun, radius_Rsun, radius_conv_core_Rsun, separation_inner_AU, mass_Msun, k_2, debug_flag);
tau_dyn_circ_yr = calculate_timescale_dynamical_tides_circ(mass_Msun, mass_conv_core_Msun, radius_Rsun, radius_conv_core_Rsun, separation_inner_AU, mass_Msun, debug_flag);
tau_tide        = calculate_tau_ZLK_over_episilon_tide(mass_Msun, mass_conv_core_Msun, radius_Rsun, radius_conv_core_Rsun, separation_inner_AU, mass_Msun, k_2, debug_flag);

if debug_flag
    separation_inner_AU
    tau_dyn_sync_yr
    tau_dyn_circ_yr
    tau_tide
end

% Create 2D meshgrid
[X, Y]  = meshgrid(m3, Pout);

% Calculate stability criteria
mock_mass                           = linspace(1,100,1001);
q_out                               = mock_mass./(mass_Msun+mass_Msun);
crit_stability_Vynatheya_a_AU       = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr   = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d    = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;
crit_stability_P_orb_d              = crit_stability_Vynatheya_P_orb_d;

% Calculate ZLK timescales
tau_ZLK(find(tau_ZLK>=10^4))=10^4;

time_tolerance_0    = 100;
idx_tau_dyn_circ    = (tau_ZLK>=tau_dyn_circ_yr-time_tolerance_0) & (tau_ZLK<=tau_dyn_circ_yr+time_tolerance_0);
time_tolerance_1    = 1;
idx_tau_dyn_sync    = (tau_ZLK>=tau_dyn_sync_yr-time_tolerance_1) & (tau_ZLK<=tau_dyn_sync_yr+time_tolerance_1);
time_tolerance_2    = 0.001;
idx_tau_tide        = (tau_ZLK>=tau_tide-time_tolerance_2) & (tau_ZLK<=tau_tide+time_tolerance_2);

t_ZLK_Liu_separations(mass_Msun, mass_Msun, 1, separation_in_AU(100./AstroConstants.yr_to_d,mass_Msun+mass_Msun+1), separation_inner_AU)

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
plot_label_png = strcat(plot_label_temp_png,'tau_ZLK.png');
plot_label_pdf = strcat(plot_label_temp_pdf,'tau_ZLK.pdf');

sz = 145.0;
fs=21;
fs2=21;
lw=2.0;
chosen_aspect = 1.5;
step_number = 15;

dark_grey   = 0.3.*[1 1 1];

xLims       = [1 100];
xLimsLabels = {'1','10','100'};
yLims       = [1 1000];
yLimsLabels = {'1','10','100','1000'};

leveler = 10;

clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')

hold on        
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink(1000)))
cbar.Label.Interpreter = 'latex';
cbar.XTick = [-4:1:4];

% ZLK timescale
extreme_val_tau = max(abs(log10(min(min(tau_ZLK)))),abs(log10(max(max(tau_ZLK)))));
cbar.Label.String = '$\log_{10} (\tau_{\rm{ZLK}}/\rm{yr})$';
mm_tau_ZLK=mesh(X, Y, log10(tau_ZLK'),'HandleVisibility','off');
mm_tau_ZLK.FaceColor = 'flat';
scatter(10000,10000,sz,-extreme_val_tau,'HandleVisibility','off')
scatter(10000,10000,sz,extreme_val_tau,'HandleVisibility','off')
empty_region = fill3([1 100 100 1],[10 10 1000 1000 ], -2.*[1 1 1 1],'k','HandleVisibility','off');    
colormap(flip(slanCM('coolwarm',step_number)))        

if annotations_flag
    text(1.25,300,leveler,'$\tau_{\rm{circ}}=1535\ \rm{yr}$','Color','w','Fontsize',fs2,'Rotation',12.5)
    text(4,62.5,leveler,'$\tau_{\rm{sync}}=23\ \rm{yr}$','Color','k','Fontsize',fs2,'Rotation',12.5)    
    text(18,7.5,leveler,'$\tau_{\rm{aps}}=0.1\ \rm{yr}$','Color','k','Fontsize',fs2,'Rotation',12.5)        
end

plot3([min(X(idx_tau_dyn_circ')), max(X(idx_tau_dyn_circ'))],[min(Y(idx_tau_dyn_circ')), max(Y(idx_tau_dyn_circ'))],[leveler+1 leveler+1],':w','LineWidth',lw)
plot3([min(X(idx_tau_dyn_sync')), max(X(idx_tau_dyn_sync'))],[min(Y(idx_tau_dyn_sync')), max(Y(idx_tau_dyn_sync'))],[leveler leveler],'--k','LineWidth',lw)
plot3([min(X(idx_tau_tide')), max(X(idx_tau_tide'))],[min(Y(idx_tau_tide')), max(Y(idx_tau_tide'))],[leveler-1 leveler-1],'k','LineWidth',lw)

dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],dark_grey);
unstableRegion.EdgeColor = 'none';
unstableRegion.FaceColor = dark_grey;
unstableRegion.HandleVisibility = 'off';
text(1.2,1.75,leveler,'Dynamically Unstable','Color','w','Fontsize',fs)

ylabel('$P_{\rm{out}}/\rm{d}$')
xlabel('$m_3/M_{\odot}$')

pbaspect([chosen_aspect 1 1])
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

function tau_circ_dyn_yr = calculate_timescale_dynamical_tides_circ(M_Msun, Mconv_Msun, R_Rsun, Rconv_Rsun, a_AU, M_comp, debug_flag)
% Based on Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
q2 = M_comp./M_Msun;
R_AU = R_Rsun.*AstroConstants.Rsun_to_AU;
E_2 = calculate_E_2(M_Msun, R_Rsun, Rconv_Rsun, debug_flag);
inv_tau_dyn_circ = (21.0./2).*sqrt((AstroConstants.G_AUMsunyr.*M_Msun./(R_AU.^3))).*q2.*((1+q2).^(11.0/6)).*E_2.*((R_AU/a_AU).^(21.0/2)); % Eq. 41 of Hurley, Pols & Tout (2002)
tau_circ_dyn_yr = 1./inv_tau_dyn_circ;
end

function tau_sync_dyn_yr = calculate_timescale_dynamical_tides_sync(M_Msun, Mconv_Msun, R_Rsun, Rconv_Rsun, a_AU, M_comp, k, debug_flag)
% Based on Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
q2 = M_comp./M_Msun;
R_AU = R_Rsun.*AstroConstants.Rsun_to_AU;
beta = 0.7.*k^0.224; % Fit following Claret & Gimenez (1989).
E_2 = calculate_E_2(M_Msun, R_Rsun, Rconv_Rsun, debug_flag);

% Eq. 44 of Hurley, Pols & Tout (2002)
inv_tau_dyn_sync = (5*(2.0.^(5.0/3))).*sqrt((AstroConstants.G_AUMsunyr.*M_Msun./(R_AU.^3))).*(1./(beta.*beta)).*q2.*q2.*((1+q2).^(5.0/6)).*E_2.*((R_AU/a_AU).^(17.0/2));
tau_sync_dyn_yr = 1./inv_tau_dyn_sync;

end

function E_2 = calculate_E_2(M_Msun, R_Rsun, Rconv_Rsun, debug_flag)

list_rad = {'Zahn','Yoon','Qin'};
[indx_rad, tf_rad] = listdlg('ListString',list_rad);

if indx_rad==1
% Based on Zahn (1975)
% As implemented in Hurley, Pols & Tout (2002)
% https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract
E_2 = 1.592*(10^-9).*(M_Msun.^(2.84)); % Eq. 43 of Hurley, Pols & Tout (2002)

elseif indx_rad==2
% Based on Yoon et al. (2010)
% https://ui.adsabs.harvard.edu/abs/2010ApJ...725..940Y/abstract
% Equation (5)
E_2 = (10.0^-1.37)*(Rconv_Rsun./R_Rsun).^8;
elseif indx_rad==3
% Based on Qin et al. (2018)
% https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..28Q/abstract
E_2 = (10.0^-0.42)*(Rconv_Rsun./R_Rsun).^7.5;
end

if debug_flag
    list_rad{indx_rad}
    E_2
end

end

function tau_tide = calculate_tau_ZLK_over_episilon_tide(M_Msun, Mconv_Msun, R_Rsun, Rconv_Rsun, a_AU, M_comp, k, debug_flag)
% tau_ZLK/epsilon_tide
Rconv_AU = R_Rsun.*AstroConstants.Rsun_to_AU;
tau_tide = (1.0/60).*(M_comp./M_Msun).*(1./sqrt(AstroConstants.G_AUMsunyr.*(M_Msun+M_comp))).*(a_AU.^(13.0/2)).*(1./k).*(1./Rconv_AU.^5.0);

end

function t_ZLK_yr = t_ZLK_Liu_separations(m1_Msun, m2_Msun, m3_Msun, aout_AU, ain_AU)

numerator   = (sqrt(m1_Msun+m2_Msun)).*(aout_AU.^3);
denominator = (sqrt(AstroConstants.G_AUMsunyr)).*m3_Msun.*(ain_AU.^(3.0/2)); 

t_ZLK_yr    = numerator./denominator;

end