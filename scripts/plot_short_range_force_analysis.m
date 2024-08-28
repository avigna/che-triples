function plot_short_range_force_analysis(debug_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% save_flag: bool

% TIC 470710327
P_out_TIC = 52.04;
m_out_TIC = 15.5;

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));
P_out_eta_unity_days    = @(m, P_d, m3) (2^-4).*m.*m.*(1./(m3.^3.0)).*(2*m+m3)*P_d;

% Stability criteria
stability_MA01 = @(a_in, q_out, e_out, i_mut) 2.8*((1+q_out).*((1+e_out)./sqrt((1-e_out))).^(2.0/5)).*(1-(0.3.*i_mut/pi));
Y_crit = @(q_out, e_in_tilde, e_out, i_mut) 2.4 * ((1 + q_out) / ((1 + e_in_tilde) * (1 - e_out)^(1/2)))^(2/5) ...
    * (((1 - 0.2 * e_in_tilde + e_out) / 8) * (cos(i_mut) - 1) + 1);

% DATA
% filename            = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA_GR_Tides.mat';
title_string        = 'm_1=m_2=55 Msun, P_{orb}=1 d, Z=0.0001';
mass_Msun           = 55;
radius_Rsun         = 7;
orbital_period_days = 1.0; 

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
list_plot = {'Test particle','ZKL','SA','SA+GR','SA+Tides','GR','GR+Tides','Tides','SA+GR+Tides'};
[indx_plot, tf_plot] = listdlg('ListString',list_plot);

if indx_plot==1
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_test_particle.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_test_particle.pdf';
    title_legend    = '$\bf ZKL\ test\ particle\ limit$'
elseif indx_plot==2
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_ZKL.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_ZKL.pdf';
    title_legend    = '$\bf ZKL$'
elseif indx_plot==3
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA.pdf';   
    title_legend    = '$\bf ZKL+SA$'    
elseif indx_plot==4
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA_GR.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_GR.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_GR.pdf';
    title_legend    = '$\bf ZKL+SA+GR$'    
elseif indx_plot==5
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_Tides.pdf';
    title_legend    = '$\bf ZKL+SA+Tides$'    
elseif indx_plot==6
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_GR.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_GR.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_GR.pdf';
    title_legend    = '$\bf ZKL+GR$'    
elseif indx_plot==7
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_GR_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_GR_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_GR_Tides.pdf';
    title_legend    = '$\bf ZKL+GR+Tides$'    
elseif indx_plot==8
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_Tides.pdf';
    title_legend    = '$\bf ZKL+Tides$'    
elseif indx_plot==9
    filename        = '../data/dynamics/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA_GR_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_GR_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_GR_Tides.pdf';
    title_legend    = '$\bf ZKL+SA+GR+Tides$'    
else
    warning("Odd choice.")
end  

M                               = load(filename);
m3                              = M.m3;
Pout                            = M.p2;
e_max                           = M.eccs;
[X, Y]                          = meshgrid(m3, Pout);

% Some quantities
num_dots = 1000;
% separation_inner_AU = separation_inner/AstroConstants.AU_to_Rsun;
Minner = 2*mass_Msun;
Pdays = [1 5 6 10 30 50 100 500 1000];
Pyears = Pdays./AstroConstants.yr_to_d;

if debug_flag
    Pdays(2)
    Pdays(3)
    Pdays(5)
end

% Calculate eta=1
mock_mass = linspace(1,100,1001);
YY = P_out_eta_unity_days(mass_Msun,orbital_period_days,mock_mass);

q_out = mock_mass./(mass_Msun+mass_Msun);
crit_stability_AU       = stability_MA01(separation_inner_AU, q_out, 0.0, 0.0).*separation_inner_AU;
crit_stability_P_orb_yr = orbital_period_yr(crit_stability_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_P_orb_d  = crit_stability_P_orb_yr.*AstroConstants.yr_to_d;

% Plot
solar=char(9737);
sz = 145.0;
lw=2.0;
fs=18;

lines1 = [         0    0.4470    0.7410];
lines2 = [    0.8500    0.3250    0.0980];
lines3 = [    0.9290    0.6940    0.1250];
lines4 = [    0.4940    0.1840    0.5560];
lines5 = [    0.4660    0.6740    0.1880];
lines6 = [    0.3010    0.7450    0.9330];
lines7 = [    0.6350    0.0780    0.1840];
dark_grey   = 0.3.*[1 1 1];
light_grey  = 0.6.*[1 1 1];    
grey   = 0.5.*[1 1 1];

pos_hor_6_d = 9.5;
pos_ver_6_d = 0.41;
pos_hor_30_d = 3.5;
pos_ver_30_d = 1.1;

pos_hor     = 1;
pos_wide    = 6.0;
pos_close   = 0.5;
pos_unst    = 2;

xLims       = [1 100];
xLimsLabels = {'1','10','100'};
yLims       = [1 1000];
yLimsLabels = {'1','10','100','1000'};

dark_grey   = 0.3.*[1 1 1];
light_grey  = 0.6.*[1 1 1];

leveler     = 2;


clf
set(gca,'defaulttextinterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')
hold on
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink))
cbar.Label.String = '$e_{\rm{max}}$';
cbar.Label.Interpreter = 'latex';

if indx_plot==1
    e_max_ZKL_test_particle = ones(size(e_max));
    mm=mesh(X, Y, e_max_ZKL_test_particle');
else
    mm=mesh(X, Y, e_max');
end

mm.FaceColor = 'flat';
scatter(1000,1000,sz,0,'HandleVisibility','off')
scatter(1000,1000,sz,1,'HandleVisibility','off')
% plot(mock_mass, crit_stability_P_orb_d,'-','LineWidth',lw,'Color',lines6)

% text(pos_hor,pos_unst,'Dynamically Unstable','Color',lines6,'Fontsize',fs)
% plot3(mock_mass, crit_stability_P_orb_d,2.*ones(size(crit_stability_P_orb_d)),'-','LineWidth',lw,'Color',lines6)
dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion_Mardling = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],grey);
unstableRegion_Mardling.EdgeColor = 'none';
unstableRegion_Mardling.FaceColor = grey;
unstableRegion_Mardling.HandleVisibility = 'off';
text(1.2,2,leveler,'Dynamically Unstable','Color','w','Fontsize',fs) 

title(title_legend)

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
ax1.XTickLabels = xLimsLabels;
ax1.YTickLabels = yLimsLabels;
ax1.TickLabelInterpreter = 'latex';
set(cbar,'TickLabelInterpreter','latex')
box on

colormap(flip(pink(100)))

% Save
if save_flag
    print(gcf,plot_label_png,'-dpng','-r250');
    print(gcf,plot_label_pdf,'-dpdf');    
end

end