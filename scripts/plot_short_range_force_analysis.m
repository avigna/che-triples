function plot_short_range_force_analysis(debug_flag, annotations_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% save_flag: bool

% Functions
% Kepler's law
separation_in_AU        = @(P_yr, M_Msun) (M_Msun.*P_yr.*P_yr).^(1.0/3);
orbital_period_yr       = @(a_AU, M_Msun) sqrt((a_AU.^3.0)./(M_Msun));

% Stability criteria
a_out_stability_Vynatheya_circular    = @(a_in, q_out, i_mut) (2.4*((1+q_out).^(2.0/5)).*(((cos(i_mut)-1)./8.0)+1)).*a_in;

% DATA
mass_Msun               = 55;
orbital_period_days     = 1.1; 

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
    title_legend    = '$\rm d)\ ZLK+SA+GR$'    
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
    title_legend    = '$\rm e)\ ZLK+DA+Tides$'    
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
mock_mass = linspace(1,100,1001);

q_out = mock_mass./(mass_Msun+mass_Msun);
crit_stability_Vynatheya_a_AU      = a_out_stability_Vynatheya_circular(separation_inner_AU, q_out, 0.0);
crit_stability_Vynatheya_P_orb_yr  = orbital_period_yr(crit_stability_Vynatheya_a_AU, mass_Msun+mass_Msun+mock_mass);
crit_stability_Vynatheya_P_orb_d   = crit_stability_Vynatheya_P_orb_yr.*AstroConstants.yr_to_d;
crit_stability_P_orb_d = crit_stability_Vynatheya_P_orb_d;

% PLOT
sz = 145.0;
fs=21;
chosen_aspect = 1.5;

% Colors
dark_grey   = 0.3.*[1 1 1];

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

% Save
if save_flag
    print(gcf,plot_label_png,'-dpng','-r250');
    print(gcf,plot_label_pdf,'-dpdf');    
end

end