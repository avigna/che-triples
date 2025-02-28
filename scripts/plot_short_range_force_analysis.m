function plot_short_range_force_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% save_flag: bool

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));
P_out_eta_unity_days    = @(m, P_d, m3) (2^-4).*m.*m.*(1./(m3.^3.0)).*(2*m+m3)*P_d;

% Stability criteria
a_out_stability_Vynatheya_circular    = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
% # 0.1 Z_SMC
% idx_0_1_Z_SMC = [47, 1784, 2132, 2504]
% age_yr_0_1_Z_SMC = [5.329669e+03, 5.120677e+06, 5.307573e+06, 5316369.558564772]
% apsidal_constant_k2_0_1_Z_SMC = [0.023533, 0.008881, 0.004897, 0.0007136587012452427]
% period_days_0_1_Z_SMC = [1.099654, 1.450562, 1.668519, 1.6824947859764587]
% radius_Rsun_0_1_Z_SMC = [7.577242, 2.401224, 1.939817, 0.8534730543749688]
% mass_conv_core_Msun_0_1_Z_SMC = [37.805430, 41.412994, 38.487511, 0.0]
% mass_Msun_0_1_Z_SMC = [54.999836, 43.976249, 43.793192, 43.79319230382309]
title_string        = 'm_1=m_2=55 Msun, P_{orb}=1.1 d, Z=0.00035';
mass_Msun               = 55;
% mass_conv_core_Msun     = 30.1;
% radius_conv_core_Rsun   = 3.2;
radius_Rsun             = round(7.577242,2);
orbital_period_days     = 1.1; 
% luminosity_Lsun         = 10^5.5; % Pending value

% Calculate extra values
orbital_period_year = orbital_period_days./AstroConstants.yr_to_d;
separation_inner_AU = separation_in_AU(orbital_period_year,mass_Msun+mass_Msun);

% Print values
if debug_flag
    fprintf('Mass_1 = Mass 2 = %f',mass_Msun)
    fprintf('\n')
    fprintf('a_{inner} = %f', separation_inner)
    fprintf('\n')
end

% Choose plot
list_plot = {'Test particle','ZLK','SA','SA+GR','SA+Tides','GR','GR+Tides','Tides','SA+GR+Tides'};
[indx_plot, tf_plot] = listdlg('ListString',list_plot);

color_title = 'k';
if indx_plot==1
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654.mat';
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_test_particle_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_test_particle_low_Z.pdf';
    title_legend    = '$\rm ZLK\ test\ particle\ limit$'
    color_title = 'w';
elseif indx_plot==2
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_ZLK_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_ZLK_low_Z.pdf';
    title_legend    = '$\rm a)\ ZLK+DA$'
    color_title = 'w';    
elseif indx_plot==3
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA.mat';
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_SA_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_SA_low_Z.pdf';   
    title_legend    = '$\rm b)\ ZLK+SA$'    
    color_title = 'w';    
elseif indx_plot==4
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR.mat';
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_SA_GR_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_SA_GR_low_Z.pdf';
    title_legend    = '$\rm e)\ ZLK+SA+GR$'    
elseif indx_plot==5
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_Tides.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_SA_Tides_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_SA_Tides_low_Z.pdf';
    title_legend    = '$\rm f)\ ZLK+SA+Tides$'    
elseif indx_plot==6
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_GR.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_GR_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_GR_low_Z.pdf';
    title_legend    = '$\rm c)\ ZLK+DA+GR$'    
elseif indx_plot==7
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_GR_Tides.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_GR_Tides_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_GR_Tides_low_Z.pdf';
    title_legend    = '$\rm g)\ ZLK+DA+GR+Tides$'    
elseif indx_plot==8
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_Tides.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_Tides_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_Tides_low_Z.pdf';
    title_legend    = '$\rm d)\ ZLK+DA+Tides$'    
elseif indx_plot==9
    filename        = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR_Tides.mat';    
    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-short-range-force-analysis_SA_GR_Tides_low_Z.png';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-short-range-force-analysis_SA_GR_Tides_low_Z.pdf';
    title_legend    = '$\rm h)\ ZLK+SA+GR+Tides$'    
else
    warning("Odd choice.")
end  

M                               = load(filename);
m3                              = M.m3;
Pout                            = M.p2;
e_max                           = M.eccs;
[X, Y]                          = meshgrid(m3, Pout);

% Some quantities
% num_dots = 1000;
% separation_inner_AU = separation_inner/AstroConstants.AU_to_Rsun;
% Minner = 2*mass_Msun;
Pdays = [1 5 6 10 30 50 100 500 1000];
% Pyears = Pdays./AstroConstants.yr_to_d;

if debug_flag
    Pdays(2)
    Pdays(3)
    Pdays(5)
end

% Calculate eta=1
mock_mass = linspace(1,100,1001);
% YY = P_out_eta_unity_days(mass_Msun,orbital_period_days,mock_mass);

q_out = mock_mass./(mass_Msun+mass_Msun);
crit_stability_Vynatheya_a_AU      = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr  = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d   = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;
crit_stability_P_orb_d = crit_stability_Vynatheya_P_orb_d;

% Filtering mergers
e_lim_val           = 1-(2*radius_Rsun/(separation_inner_AU.*AstroConstants.AU_to_Rsun));
e_tol               = 0.0075;
idx_of_eccentricity = find((e_max'<=e_lim_val+e_tol)&(e_max'>=e_lim_val-e_tol));

% PLOT
sz = 145.0;
fs=21;
chosen_aspect = 1.5;

% Colors
lines1 = [0    0.4470    0.7410];

dark_grey   = 0.3.*[1 1 1];
light_grey  = 0.6.*[1 1 1];    
grey   = 0.5.*[1 1 1];

% pos_hor_6_d = 9.5;
% pos_ver_6_d = 0.41;
% pos_hor_30_d = 3.5;
% pos_ver_30_d = 1.1;
% 
% pos_hor     = 1;
% pos_wide    = 6.0;
% pos_close   = 0.5;
% pos_unst    = 2;

xLims       = [1 100];
xLimsLabels = {'1','10','100'};
yLims       = [1 1000];
yLimsLabels = {'1','10','100','1000'};

text_x_coord = 1.2;

leveler     = 2;

clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')
hold on
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink(100)))
cbar.Label.String = '$e_{\rm{max}}$';
cbar.Label.Interpreter = 'latex';
cbar.XTick = [0:0.2:1];

if indx_plot==1
    e_max_ZLK_test_particle = ones(size(e_max));
    mm=mesh(X, Y, e_max_ZLK_test_particle');
else
    mm=mesh(X, Y, e_max');
end

mm.FaceColor = 'flat';
scatter(1000,1000,sz,0,'HandleVisibility','off')
scatter(1000,1000,sz,1,'HandleVisibility','off')

dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],dark_grey);
unstableRegion.EdgeColor = 'none';
unstableRegion.FaceColor = dark_grey;
unstableRegion.HandleVisibility = 'off';
text(text_x_coord,2,leveler+1,'Dynamically Unstable','Color','w','Fontsize',fs)

text(1.2,500,leveler, title_legend,'Color',color_title,'Fontsize',fs)

if annotations_flag
 % plot3(mock_mass,YY,ones(size(YY)),'Color',lines1,'Linewidth',lw)
    % scatter3(X(idx_of_eccentricity),Y(idx_of_eccentricity),1000.*ones(size(X(idx_of_eccentricity))),20,lines1,'s','filled')
end

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
ax1.XTickLabels = xLimsLabels;
ax1.YTickLabels = yLimsLabels;
ax1.TickLabelInterpreter = 'latex';
set(cbar,'TickLabelInterpreter','latex')


box on

% colormap(flip(pink(100)))

% Save
if save_flag
    print(gcf,plot_label_png,'-dpng','-r250');
    print(gcf,plot_label_pdf,'-dpdf');    
end

end