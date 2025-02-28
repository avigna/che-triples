function plot_close_triple_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% annotations_flag: bool
% save_flag: bool

% TIC 470710327
P_out_TIC = 52.04;
m_out_TIC = 15.5;

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));
% P_out_eta_unity_days    = @(m, P_in_d, m3) (2^-4).*m.*m.*(1./(m3.^3.0)).*(2*m+m3)*P_d;

% Stability criteria
a_out_stability_Vynatheya_circular    = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
% # Z_SMC
% idx_Z_SMC = [39, 184, 1340, 2286]
% age_yr_Z_SMC = [8.392999e+03, 2.697549e+06, 5.310639e+06, 5.683500e+06]
% apsidal_constant_k2_Z_SMC = [0.018226, 0.011376, 0.006834, 0.001726]
% period_days_Z_SMC = [1.100177, 1.087004, 2.220218, 4.525685]
% radius_Rsun_Z_SMC = [9.038323, 10.171930, 2.016422, 0.448425]
% mass_conv_core_Z_SMC = [37.751202, 38.110672, 19.155549, 0.000000]
% mass_Msun_Z_SMC = [54.997851, 51.413071, 26.583868, 18.619564]

% # 0.1 Z_SMC
% idx_0_1_Z_SMC = [47, 1784, 2132, 2504]
% age_yr_0_1_Z_SMC = [5.329669e+03, 5.120677e+06, 5.307573e+06, 5316369.558564772]
% apsidal_constant_k2_0_1_Z_SMC = [0.023533, 0.008881, 0.004897, 0.0007136587012452427]
% period_days_0_1_Z_SMC = [1.099654, 1.450562, 1.668519, 1.6824947859764587]
% radius_Rsun_0_1_Z_SMC = [7.577242, 2.401224, 1.939817, 0.8534730543749688]
% mass_conv_core_Msun_0_1_Z_SMC = [37.805430, 41.412994, 38.487511, 0.0]
% mass_Msun_0_1_Z_SMC = [54.999836, 43.976249, 43.793192, 43.79319230382309]
% Choose
list_binary = {'55+55 CHE ZAMS low Z','55+55 CHE HeMS low Z','55+55 CHE end low Z','55+55 CHE ZAMS high Z','55+55 CHE HeMS high Z','55+55 CHE end high Z','TIC 470710327 full','TIC 470710327 GR'};
[indx_binary, tf_binary] = listdlg('ListString',list_binary);


if indx_binary==1
    filename            = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR_Tides.mat';
    title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.00035/ \rm{at\ ZAMS}$'];
    plot_label_temp_png = '../plots/png/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_low_Z_';
    plot_label_temp_pdf = '../plots/pdf/55_Msun_low_Z/CHE-ZAMS-close-triple-analysis_low_Z_';    
    mass_Msun               = 55;
    radius_Rsun             = round(7.577242,2);
    orbital_period_days     = 1.1; 
elseif indx_binary==2
    filename            = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=43.976249_Porb=1.450562_SA_GR_Tides.mat';
    title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.00035/ \rm{at\ HeMS}$'];
    plot_label_temp_png = '../plots/png/55_Msun_low_Z/CHE-HeMS-close-triple-analysis_low_Z_';
    plot_label_temp_pdf = '../plots/pdf/55_Msun_low_Z/CHE-HeMS-close-triple-analysis_low_Z_';    
    mass_Msun               = round(43.976249,2)
    radius_Rsun             = round(2.401224,2)
    orbital_period_days     = round(1.450562,1); 
elseif indx_binary==3
    filename            = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=43.79319230382309_Porb=1.6824947859764587_SA_GR.mat';
    title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.00035/ \rm{at\ end}$'];
    plot_label_temp_png = '../plots/png/55_Msun_low_Z/CHE-end-close-triple-analysis_low_Z_';
    plot_label_temp_pdf = '../plots/pdf/55_Msun_low_Z/CHE-end-close-triple-analysis_low_Z_';    
    mass_Msun               = round(43.79319230382309,2)
    orbital_period_days     = round(1.6824947859764587,2);     
elseif indx_binary==4
    filename            = '../data/dynamics/55_Msun_Z_SMC/triple_Z=0.0035_CHE=1_M1=M2=54.997851_Porb=1.100177_SA_GR_Tides.mat';
    title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.0035/ \rm{at\ HeMS}$'];
    plot_label_temp_png = '../plots/png/55_Msun_high_Z/CHE-ZAMS-close-triple-analysis_high_Z_';
    plot_label_temp_pdf = '../plots/pdf/55_Msun_high_Z/CHE-ZAMS-close-triple-analysis_high_Z_';    
    mass_Msun               = round(54.997851,2)
    radius_Rsun             = round(9.038323,2)
    orbital_period_days     = round(1.100177,1); 
elseif indx_binary==5
    filename            = '../data/dynamics/55_Msun_Z_SMC/triple_Z=0.0035_CHE=1_M1=M2=26.583868_Porb=2.220218_SA_GR_Tides.mat';
    title_string        = ['$m_1=m_2=55\ M_{\odot},\ P_{\rm{orb}}=1.1\ d,\ Z=0.0035/ \rm{at\ HeMS}$'];
    plot_label_temp_png = '../plots/png/55_Msun_high_Z/CHE-HeMS-close-triple-analysis_high_Z_';
    plot_label_temp_pdf = '../plots/pdf/55_Msun_high_Z/CHE-HeMS-close-triple-analysis_high_Z_';    
    mass_Msun               = round(26.583868,2)
    radius_Rsun             = round(2.016422,2)
    orbital_period_days     = round(2.220218,1); 
elseif indx_binary==6
     
elseif indx_binary==7
    % filename            = '../data/dynamics/6_Msun/triple_Z=0.142_CHE=0_M1=M2=6_Porb=1.1_SA_GR_Tides.mat';
    % title_string        = '$m_1=m_2=6\ M_{\odot}, P_{\rm{orb}}=1.1\ d,\ Z=0.0142$';
    % plot_label_temp_png = '../plots/png/TIC-full-close-triple-analysis_';
    % plot_label_temp_pdf = '../plots/pdf/TIC-full-close-triple-analysis_';        
    % mass_Msun           = 6;
    % radius_Rsun         = 2.8;
    % orbital_period_days = 1.1;
elseif indx_binary==8
    % filename            = '../data/dynamics/6_Msun/triple_Z=0.142_CHE=0_M1=M2=6_Porb=1.1_SA_GR.mat';
    % title_string        = '$m_1=m_2=6\ M_{\odot}, P_{\rm{orb}}=1.1\ d,\ Z=0.0142$';
    % plot_label_temp_png = '../plots/png/TIC-GR-close-triple-analysis_';
    % plot_label_temp_pdf = '../plots/pdf/TIC-GR-close-triple-analysis_';        
    % mass_Msun           = 6;
    % radius_Rsun         = 2.8;
    % orbital_period_days = 1.1;    
else
    warning("Odd choice.")
end

% Load
if debug_flag==true
    M                               = load(filename)
else
    M                               = load(filename);
end

cos_inc     = M.cos_inc;
e_max       = M.eccs;
eps_SA      = M.eps_SA;
eps_GR      = M.eps_GR;
eps_Tide    = M.eps_Tide;
eta         = M.eta;
f_merger    = M.f_merger;
m3          = M.m3;
Pout        = M.p2;
tau_sec     = M.tau_sec.*AstroConstants.s_to_yr;

max(max(eps_SA))
% Calculate extra values
orbital_period_year = orbital_period_days./AstroConstants.yr_to_d;
separation_inner_AU = separation_in_AU(orbital_period_year,mass_Msun+mass_Msun);

% Create 2D meshgrid
[X, Y]                          = meshgrid(m3, Pout);

% Calculate stability
mock_mass   = linspace(1,100,1001);
q_out       = mock_mass./(mass_Msun+mass_Msun);

crit_stability_Vynatheya_a_AU       = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr   = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d    = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;
crit_stability_P_orb_d              = crit_stability_Vynatheya_P_orb_d;

% Filtering mergers
if indx_binary==3 

else
e_lim_val = 1-(2*radius_Rsun/(separation_inner_AU.*AstroConstants.AU_to_Rsun));
e_lim_tol = 0.01;
idx_of_eccentricity = find((e_max'>=e_lim_val-e_lim_tol)&(e_max'<=e_lim_val+e_lim_tol)&(Pout'>=5));

min_eccentricity = 0.001;
index_of_non_induced_eccentricity = find(e_max<min_eccentricity);
cos_inc(index_of_non_induced_eccentricity) = nan;
e_max(index_of_non_induced_eccentricity) = nan;
end

min_val_colorbar    = log10(min(min(min(min(eps_SA)),min(min(eps_GR))),min(min(eps_Tide))));
max_val_colorbar    = log10(max(max(max(max(eps_SA)),max(max(eps_GR))),max(max(eps_Tide))));
extreme_val_colorbar= max(abs(min_val_colorbar),abs(max_val_colorbar));

% Print values
if debug_flag==true
    title_string
    display(list_binary(indx_binary))
    fprintf('Mass_1 = Mass 2 = %f Msun',mass_Msun)
    fprintf('\n')
    fprintf('a_{inner} = %f au', separation_inner_AU)
    fprintf('\n')
    fprintf('P_{orb} = %f d', orbital_period_days)
    fprintf('\n')    
    fprintf('e_{lim} = %f', e_lim_val)
    fprintf('\n')    
end

% PLOT
list_plot = {'cos(i)|_{e_{max}}','e_max','\eps_{GR}','\eps_{SA}','\eps_{Tide}','\eta','f_merger','\tau_{sec}'};
[indx_plot, tf_plot] = listdlg('ListString',list_plot);

if indx_plot==1
    plot_label_png = strcat(plot_label_temp_png,'inclination.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'inclination.pdf');
elseif indx_plot==2
    plot_label_png = strcat(plot_label_temp_png,'max_eccentricity.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'max_eccentricity.pdf');
elseif indx_plot==3
    plot_label_png = strcat(plot_label_temp_png,'eps_GR.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'eps_GR.pdf');    
elseif indx_plot==4
    plot_label_png = strcat(plot_label_temp_png,'eps_SA.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'eps_SA.pdf');
elseif indx_plot==5
    plot_label_png = strcat(plot_label_temp_png,'eps_Tides.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'eps_Tides.pdf');
elseif indx_plot==6
    plot_label_png = strcat(plot_label_temp_png,'eta.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'eta.pdf');    
elseif indx_plot==7
    plot_label_png = strcat(plot_label_temp_png,'f_merger.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'f_merger.pdf');    
elseif indx_plot==8
    plot_label_png = strcat(plot_label_temp_png,'tau_sec.png');
    plot_label_pdf = strcat(plot_label_temp_pdf,'tau_sec.pdf');
else
    warning("Odd choice.")
end  


solar=char(9737);
sz = 145.0;
sz2=25;
lw=2.0;
fs=21;
chosen_aspect = 1.5;
step_number = 15;

lines6      = [    0.3010    0.7450    0.9330];
dark_grey   = 0.3.*[1 1 1];
light_grey  = 0.6.*[1 1 1];    
grey   = 0.5.*[1 1 1];

xLims       = [1 100];
xLimsLabels = {'1','10','100'};
yLims       = [1 1000];
yLimsLabels = {'1','10','100','1000'};

extreme_val_eta = max(abs(log10(min(min(eta)))),abs(log10(max(max(eta)))));
extreme_val_tau = max(abs(log10(min(min(tau_sec)))),abs(log10(max(max(tau_sec)))));

leveler = 100;
text_x_coord = 1.2;
eps_ticks_color_bar = [-4:2:4];

clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')

hold on        
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink(1000)))
cbar.Label.Interpreter = 'latex';

if indx_plot == 1
    % inclination for maximum eccentricity
    cbar.Label.String = '$\cos(i_0)|_{e_{\rm{max}}}$';
    mm_inclination=mesh(X, Y, cos_inc','HandleVisibility','off');
    mm_inclination.FaceColor = 'flat';
    scatter(10000,10000,sz,-1,'HandleVisibility','off')
    scatter(10000,10000,sz,0,'HandleVisibility','off') 
    empty_region = fill3([1 100 100 1],[10 10 1000 1000 ], -2.*[1 1 1 1],'k');
    % colormap(flip(slanCM('magma',10)))
    % colormap(flip(pink))
    cbar.XTick = [-1:0.2:0];
elseif indx_plot == 2
    % maximum eccentricity
    cbar.Label.String = '$e_{\rm{max}}$';
    mm_eccentricity=mesh(X, Y, e_max','HandleVisibility','off');
    mm_eccentricity.FaceColor = 'flat';
    scatter(10000,10000,sz,0,'HandleVisibility','off')
    scatter(10000,10000,sz,1,'HandleVisibility','off')
    empty_region = fill3([1 100 100 1],[10 10 1000 1000 ], min_eccentricity.*[1 1 1 1],'k');
    % colormap(flip(slanCM('magma',10)))
    % colormap(flip(pink))    
elseif indx_plot == 3
    % epsilon GR
    cbar.Label.String = '$\log_{10} \epsilon_{\rm{GR}}$';
    mm_GR=mesh(X, Y, log10(eps_GR'));
    mm_GR.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_colorbar,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_colorbar,'HandleVisibility','off')
    % colormap(flip(slanCM('vik')))
    colormap(flip(slanCM('vik',step_number)))        
    cbar.XTick = eps_ticks_color_bar;
elseif indx_plot == 4
    % epsilon SA
    cbar.Label.String = '$\log_{10} \epsilon_{\rm{SA}}$';
    mm_SA=mesh(X, Y, log10(eps_SA'));
    mm_SA.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_colorbar,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_colorbar,'HandleVisibility','off')
    % colormap(flip(slanCM('vik')))
    colormap(flip(slanCM('vik',step_number)))    
    cbar.XTick = eps_ticks_color_bar;
elseif indx_plot == 5
    % epsilon Tides
    cbar.Label.String = '$\log_{10} \epsilon_{\rm{Tide}}$';
    mm_Tide=mesh(X, Y, log10(eps_Tide'));
    mm_Tide.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_colorbar,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_colorbar,'HandleVisibility','off')
    % colormap(flip(slanCM('vik')))
    colormap(flip(slanCM('vik',step_number)))        
    cbar.XTick = eps_ticks_color_bar;    
elseif indx_plot == 6
    % eta
    % cbar.Label.String = '$\log_{10} \eta := L_{\rm{in}}/L_{\rm{out}}$';
    cbar.Label.String = '$\log_{10} \eta$';    
    mm_eta=mesh(X, Y, log10(eta'));
    mm_eta.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_eta,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_eta,'HandleVisibility','off')
    colormap(flip(slanCM('vik')))    
    % colormap(flip(slanCM('vik',step_number)))       
elseif indx_plot == 7
    % f_merger
    cbar.Label.String = '$f_{\rm{merger}}$';
    mm_f_merger=mesh(X, Y, f_merger','HandleVisibility','off');
    mm_f_merger.FaceColor = 'flat';
    scatter(10000,10000,sz,0,'HandleVisibility','off')
    scatter(10000,10000,sz,1,'HandleVisibility','off')    
    cbar.XTick = [0:0.2:1];    
    % colormap(flip(pink(10)))
elseif indx_plot == 8    
    % secular timescale
    cbar.Label.String = '$\log_{10} (\tau_{\rm{sec}}/\rm{yr})$';
    mm_tau_sec=mesh(X, Y, log10(tau_sec'));
    mm_tau_sec.FaceColor = 'flat';
    scatter(10000,10000,sz,-extreme_val_tau,'HandleVisibility','off')
    scatter(10000,10000,sz,extreme_val_tau,'HandleVisibility','off')
    colormap(flip(slanCM('vik')))        
    colormap(flip(slanCM('vik')))
    cbar.XTick = [-4:2:4];
else    
    warning("Odd choice.")    
end

if annotations_flag
    % title(title_string)

    empty_region.EdgeColor = 'none';
    empty_region.FaceColor = light_grey;
    empty_region.HandleVisibility = 'off';  

    text(16,500,1,'Damped ZLK','Color','w','Fontsize',fs)            

    if indx_binary==1 
        if indx_plot==1
            text(1.2,500,1,'ZAMS','Color','k','Fontsize',fs)

        end
        if indx_plot==7
            text(1.2,500,1,'ZAMS','Color','k','Fontsize',fs)
        end
    end

    if indx_binary==2 
        if indx_plot==1
            text(1.2,500,1,'5.12 Myr','Color','k','Fontsize',fs)
        end
        if indx_plot==7
            text(1.2,500,1,'5.12 Myr','Color','k','Fontsize',fs)
        end
    end    
        
    if indx_binary==7 | indx_binary==8
        plot3(m_out_TIC,P_out_TIC,1,'or','MarkerSize',10,'MarkerFaceColor','r')
        text(m_out_TIC+2,P_out_TIC-4,1,'TIC 470710327','Color','w','Fontsize',fs)
    end


end


dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],dark_grey);
unstableRegion.EdgeColor = 'none';
unstableRegion.FaceColor = dark_grey;
unstableRegion.HandleVisibility = 'off';
text(text_x_coord,2,leveler+1,'Dynamically Unstable','Color','w','Fontsize',fs)

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