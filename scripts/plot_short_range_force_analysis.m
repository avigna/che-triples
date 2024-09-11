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
% a_out_stability_MA01 = @(a_in, q_out, e_out, i_mut) (2.8*((1+q_out).*((1+e_out)./sqrt((1-e_out))).^(2.0/5)).*(1-(0.3.*i_mut/pi))).*a_in;
% a_out_stability_Vynatheya_simple = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;
% Y_crit = @(q_out, e_in_tilde, e_out, i_mut) 2.4 * ((1 + q_out) / ((1 + e_in_tilde) * (1 - e_out)^(1/2)))^(2/5) ...
%     * (((1 - 0.2 * e_in_tilde + e_out) / 8) * (cos(i_mut) - 1) + 1);
a_out_stability_MA01_circular         = @(a_in, q_out, i_mut_rad) (2.8*((1+q_out).^(2.0/5))).*a_in;
a_out_stability_Vynatheya_circular    = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
% filename            = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=55_Porb=1_SA_GR_Tides.mat';
title_string        = 'm_1=m_2=45 Msun, P_{orb}=1 d, Z=0.0001';
mass_Msun           = 45;
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
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_test_particle.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_test_particle.pdf';
    title_legend    = '$\bf ZKL\ test\ particle\ limit$'
elseif indx_plot==2
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_ZKL.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_ZKL.pdf';
    title_legend    = '$\bf a)\ ZKL$'
elseif indx_plot==3
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_SA.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA.pdf';   
    title_legend    = '$\bf b)\ ZKL+SA$'    
elseif indx_plot==4
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_SA_GR.mat';
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_GR.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_GR.pdf';
    title_legend    = '$\bf e)\ ZKL+SA+GR$'    
elseif indx_plot==5
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_SA_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_Tides.pdf';
    title_legend    = '$\bf f)\ ZKL+SA+Tides$'    
elseif indx_plot==6
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_GR.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_GR.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_GR.pdf';
    title_legend    = '$\bf c)\ ZKL+GR$'    
elseif indx_plot==7
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_GR_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_GR_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_GR_Tides.pdf';
    title_legend    = '$\bf g)\ ZKL+GR+Tides$'    
elseif indx_plot==8
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_Tides.pdf';
    title_legend    = '$\bf d)\ ZKL+Tides$'    
elseif indx_plot==9
    filename        = '../data/dynamics/45_Msun/triple_Z=0.0001_CHE=1_M1=M2=45_Porb=1_SA_GR_Tides.mat';    
    plot_label_png  = '../plots/png/CHE-short-range-force-analysis_SA_GR_Tides.png';
    plot_label_pdf  = '../plots/pdf/CHE-short-range-force-analysis_SA_GR_Tides.pdf';
    title_legend    = '$\bf h)\ ZKL+SA+GR+Tides$'    
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
crit_stability_MA01_a_AU            = a_out_stability_MA01_circular(separation_inner_AU, q_out, 0.0);
crit_stability_MA01_P_orb_yr        = orbital_period_yr(crit_stability_MA01_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_MA01_P_orb_d         = crit_stability_MA01_P_orb_yr.*AstroConstants.yr_to_d;

crit_stability_Vynatheya_a_AU      = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr  = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d   = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;

% PLOT
list_stability = {'MA01','Vynatheya'};
[indx_stability, tf_stability] = listdlg('ListString',list_stability);

if indx_stability==1
    crit_stability_P_orb_d = crit_stability_MA01_P_orb_d;
elseif indx_stability==2
    crit_stability_P_orb_d = crit_stability_Vynatheya_P_orb_d;
else
    warning("Odd choice.")
end  


solar=char(9737);
sz = 145.0;
lw=2.0;
fs=18;

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

dummy_ones = ones(size(crit_stability_P_orb_d));
unstableRegion = fill3([mock_mass fliplr(mock_mass)],[dummy_ones fliplr(crit_stability_P_orb_d)],leveler.*[dummy_ones dummy_ones],dark_grey);
unstableRegion.EdgeColor = 'none';
unstableRegion.FaceColor = dark_grey;
unstableRegion.HandleVisibility = 'off';
text(1.2,1.75,leveler,'Dynamically Unstable','Color','w','Fontsize',fs)

% title(title_legend)
text(1.2,500,leveler, title_legend,'Color','k','Fontsize',fs)

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