function plot_GW_merger_timescale_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% annotations_flag: bool
% save_flag: bool

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));

% Stability criteria
a_out_stability_Vynatheya_circular   = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
% Choose metallicity
list_plot = {'Z =0.0042','Z=0.00042'};
[indx_plot, tf_plot] = listdlg('ListString',list_plot);


if indx_plot==1
    filename_ZAMS           = '../data/dynamics/55_Msun_Z_SMC/triple_Z=0.0035_CHE=1_M1=M2=54.997851_Porb=1.100177_SA_GR_Tides.mat';    
    mass_Msun_ZAMS          = round(54.997851,2);
    radius_Rsun_ZAMS        = round(9.038323,2);
    orbital_period_days_ZAMS= round(1.100177,2); 

    filename_HeMS           = '../data/dynamics/55_Msun_Z_SMC/triple_Z=0.0035_CHE=1_M1=M2=26.583868_Porb=2.220218_SA_GR_Tides.mat';    
    mass_Msun_HeMS          = round(26.583868,2);
    radius_Rsun_HeMS        = round(2.016422,2);
    orbital_period_days_HeMS= round(2.220218,2);

    filename_end            = '../data/dynamics/55_Msun_Z_SMC/triple_Z=0.0035_CHE=1_M1=M2=18.619564_Porb=4.525685_SA_GR.mat';
    mass_Msun_end           = round(18.619564,2);
    orbital_period_days_end = round(4.525685,2);    

    plot_label_png          = '../plots/png/55_Msun_high_Z/CHE-GW-timescale-analysis_high_Z';
    plot_label_pdf          = '../plots/pdf/55_Msun_high_Z/CHE-GW-timescale-analysis_high_Z';    
elseif indx_plot==2
    filename_ZAMS   = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=54.999836_Porb=1.099654_SA_GR_Tides.mat';
    mass_Msun_ZAMS          = round(54.999836,2);
    radius_Rsun_ZAMS        = round(7.577242,2);
    orbital_period_days_ZAMS= round(1.099654,2); 

    filename_HeMS   = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=43.976249_Porb=1.450562_SA_GR_Tides.mat';
    mass_Msun_HeMS          = round(43.976249,2);
    radius_Rsun_HeMS        = round(2.401224,2);
    orbital_period_days_HeMS= round(1.450562,2);

    filename_end    = '../data/dynamics/55_Msun_low_Z/triple_Z=0.00035_CHE=1_M1=M2=43.79319230382309_Porb=1.6824947859764587_SA_GR'; 
    
    mass_Msun_end           = round(43.79319230382309,2);
    orbital_period_days_end = round(1.6824947859764587,2);    

    plot_label_png  = '../plots/png/55_Msun_low_Z/CHE-GW-timescale-analysis_low_Z';
    plot_label_pdf  = '../plots/pdf/55_Msun_low_Z/CHE-GW-timescale-analysis_low_Z';    
else
    warning("Odd choice.")
end  
    
% Load '55+55 CHE' at ZAMS
M_ZAMS                      = load(filename_ZAMS);
e_max_ZAMS                  = M_ZAMS.eccs;
m3_ZAMS                     = M_ZAMS.m3;
Pout_ZAMS                   = M_ZAMS.p2;

% Load '55+55 CHE' at HeMS
M_HeMS                      = load(filename_HeMS);
e_max_HeMS                  = M_HeMS.eccs;
m3_HeMS                     = M_HeMS.m3;
Pout_HeMS                   = M_HeMS.p2;

% Load '55+55 CHE' at end
M_end_GR                    = load(filename_end);
e_max_end_GR                = M_end_GR.eccs;
m3_end_GR                   = M_end_GR.m3;
Pout_end_GR                 = M_end_GR.p2;

% Calculate extra values
orbital_period_year_ZAMS    = orbital_period_days_ZAMS./AstroConstants.yr_to_d;
separation_inner_AU_ZAMS    = separation_in_AU(orbital_period_year_ZAMS,mass_Msun_ZAMS+mass_Msun_ZAMS);

orbital_period_year_HeMS    = orbital_period_days_HeMS./AstroConstants.yr_to_d;
separation_inner_AU_HeMS    = separation_in_AU(orbital_period_year_HeMS,mass_Msun_HeMS+mass_Msun_HeMS);

orbital_period_year_end     = orbital_period_days_end./AstroConstants.yr_to_d;
separation_inner_AU_end     = separation_in_AU(orbital_period_year_end,mass_Msun_end+mass_Msun_end);

T_c_s_end                   = calculate_merger_time_circular_binary(separation_inner_AU_end,mass_Msun_end,mass_Msun_end);
T_c_yr_end                  = T_c_s_end.*AstroConstants.s_to_yr;

% Based on eq. 48 of Liu & Lai (2018)
% https://ui.adsabs.harvard.edu/abs/2018ApJ...863...68L/abstract
tau_ZKL_GW_yr_end           = T_c_yr_end.*((1-e_max_end_GR.^2).^3.0);

% Find less crazy values
tau_ZKL_GW_yr_end(find(tau_ZKL_GW_yr_end<=1)) = 1;
tau_ZKL_GW_yr_end(find(tau_ZKL_GW_yr_end>=10^10)) = 10^10;


% Create 2D meshgrid
[X_ZAMS, Y_ZAMS]            = meshgrid(m3_ZAMS, Pout_ZAMS);
[X_HeMS, Y_HeMS]            = meshgrid(m3_HeMS, Pout_HeMS);
[X_end_GR, Y_end_GR]        = meshgrid(m3_end_GR, Pout_end_GR);

% Calculate stability
mock_mass                   = linspace(1,100,1001);
q_out_ZAMS                  = mock_mass./(mass_Msun_ZAMS+mass_Msun_ZAMS);
q_out_HeMS                  = mock_mass./(mass_Msun_HeMS+mass_Msun_HeMS);
q_out_end                   = mock_mass./(mass_Msun_end+mass_Msun_end);

% Stability at ZAMS
crit_stability_Vynatheya_a_AU_ZAMS       = a_out_stability_Vynatheya_circular(separation_inner_AU_ZAMS, q_out_ZAMS, 0.0);
crit_stability_Vynatheya_P_orb_yr_ZAMS   = orbital_period_yr(crit_stability_Vynatheya_a_AU_ZAMS, mass_Msun_ZAMS+mass_Msun_ZAMS+mock_mass);
crit_stability_Vynatheya_P_orb_d_ZAMS    = crit_stability_Vynatheya_P_orb_yr_ZAMS.*AstroConstants.yr_to_d;

% Stability at HeMS
crit_stability_Vynatheya_a_AU_HeMS        = a_out_stability_Vynatheya_circular(separation_inner_AU_HeMS, q_out_HeMS, 0.0);
crit_stability_Vynatheya_P_orb_yr_HeMS    = orbital_period_yr(crit_stability_Vynatheya_a_AU_HeMS, mass_Msun_HeMS+mass_Msun_HeMS+mock_mass);
crit_stability_Vynatheya_P_orb_d_HeMS     = crit_stability_Vynatheya_P_orb_yr_HeMS.*AstroConstants.yr_to_d;

% Filtering mergers
e_lim_val_ZAMS           = 1-(2*radius_Rsun_ZAMS/(separation_inner_AU_ZAMS.*AstroConstants.AU_to_Rsun));
idx_of_eccentricity_ZAMS = find((e_max_ZAMS'>=e_lim_val_ZAMS)&(m3_ZAMS<99)&(m3_ZAMS>1));

e_lim_val_HeMS            = 1-(2*radius_Rsun_HeMS/(separation_inner_AU_HeMS.*AstroConstants.AU_to_Rsun));
idx_of_eccentricity_HeMS  = find((e_max_HeMS'>=e_lim_val_HeMS)&(m3_HeMS<99));


% PLOT
fs=18;
sz2=20;
leveler = 11;
lw = 2.5;

% Colors
light_grey   = 0.8.*[1 1 1];
mid_grey   = 0.5.*[1 1 1];
dark_grey   = 0.2.*[1 1 1];

chosen_color_H_mergers = mid_grey;
chosen_color_He_mergers = light_grey;

text_x_coord = 1.2;

chosen_aspect   = 1.5;
step_number     = 11;

% Limits
xLims           = [1 100];
xLimsLabels     = {'1','10','100'};
yLims           = [1 3000];
yTicksValues    = [1 10 100 1000];
yLimsLabels     = {'1','10','100','1000'};

clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')

hold on        
cbar = colorbar;
cbar.Label.Interpreter = 'latex';
cbar.XTick = [0:2:10];

% GW timescale
cbar.Label.String = '$\log_{10} (\tau_{\rm{merge,min}}/\rm{yr})$';
mm_tau_sec=mesh(X_end_GR, Y_end_GR, log10(tau_ZKL_GW_yr_end'),'HandleVisibility','off');
scatter(60, 6000, 1, 10,'HandleVisibility','off');

mm_tau_sec.FaceColor = 'flat';
cbar.FontSize = fs;
colormap(slanCM('amethyst',step_number))   


scatter3(X_HeMS(idx_of_eccentricity_HeMS),Y_HeMS(idx_of_eccentricity_HeMS),leveler-1.*ones(size(X_HeMS(idx_of_eccentricity_HeMS))),sz2,chosen_color_He_mergers,'s','filled','HandleVisibility','off')
scatter3(X_ZAMS(idx_of_eccentricity_ZAMS),Y_ZAMS(idx_of_eccentricity_ZAMS),leveler.*ones(size(X_ZAMS(idx_of_eccentricity_ZAMS))),sz2,chosen_color_H_mergers,'s','filled','HandleVisibility','off')

dummy_ones = ones(size(crit_stability_Vynatheya_P_orb_d_ZAMS));
unstableRegion_ZAMS = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_Vynatheya_P_orb_d_ZAMS)],leveler+1.*[dummy_ones dummy_ones],dark_grey);
unstableRegion_ZAMS.EdgeColor = 'none';
unstableRegion_ZAMS.FaceColor = dark_grey;
unstableRegion_ZAMS.HandleVisibility = 'off';

text(text_x_coord,2,leveler+1,'Dynamically Unstable','Color','w','Fontsize',fs)
if indx_plot==1
    text(2,34,leveler,'He Merger','Color','k','Fontsize',fs)
    text(text_x_coord,8,leveler,'Main-Sequence Merger','Color','w','Fontsize',fs)
elseif indx_plot==2
    text(2.9,22,leveler,'He Merger','Color','k','Fontsize',fs)
    text(text_x_coord,7,leveler,'Main-Sequence Merger','Color','w','Fontsize',fs)
else
    warning("Odd choice.")
end  

xline(mass_Msun_ZAMS,'--','Color','r','LineWidth',lw)


legend('$m_3=m_{1,\rm{ZAMS}}$','Location','NorthWest','interpreter','latex','Box','Off')
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
ax1.YTick = yTicksValues;
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