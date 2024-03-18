function plot_triple_analysis(string_filename, debug_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% debug_flag: bool
% save_flag: bool

% Functions
aout_crit_Mardling = @(ain,Minner,M3)  2.8.*ain.*(1+(M3./Minner)).^(2.0/5);

% From Kepler's law
isoMass = @(a_out,P_out,M_bin) ((a_out.^3)*(P_out.^-2))-M_bin;
orbital_period_yr   = @(a_AU,mtotal_Msun) sqrt((a_AU.^3)./(mtotal_Msun));
separation_au       = @(P_d, mtotal_Msun) (mtotal_Msun.*(P_d./AstroConstants.yr_to_d).^2).^(1.0/3);

% DATA
filename        = strcat('../data/dynamics/',string_filename,'.mat');
title_string    = string_filename;
plot_label      = strcat('../plots/triple_analysis_',string_filename,'.png');  

% Load
M           = load(filename);
aout        = M.a2;
cos_inc     = M.cos_inc;
e_max       = M.eccs;
eps_SA      = M.eps_SA;
eps_GR      = M.eps_GR;
eps_Tides   = M.eps_Tides;
ETA         = M.ETA;
fraction    = M.fraction;
m3          = M.m3;
rps         = M.rps.*AstroConstants.AU_to_Rsun;
tau_ZKL     = M.tau_ZKL.*AstroConstants.s_to_yr;

if debug_flag
    size(m3)
    size(aout)
    size(fraction)
    size(e_max)
    size(tau_ZKL)
end

mass_Msun           = 55;           % [mass_Msun]=Msun
radius_Rsun         = 7;
Minner              = 2*mass_Msun;  % [mass_Msun]=Msun
orbital_period_inner_days = 2.0;          % []
separation_inner_AU = separation_au(orbital_period_inner_days, Minner);
orbital_period_outer_days = orbital_period_yr(aout,m3+Minner).*AstroConstants.yr_to_d;

periastron_min_Rsun = separation_inner_AU.*(1-e_max).*AstroConstants.AU_to_Rsun;

% Print values
if debug_flag
    filename
    fprintf('Mass_1 = Mass 2 = %f',mass_Msun)
    fprintf('\n')
    fprintf('a_{inner} = %f', separation_inner)
    fprintf('\n')
    fprintf('min(f_merge) = %f',min(min(fraction)))
    fprintf('\n')
    fprintf('max(f_merge) = %f',max(max(fraction)))
    fprintf('\n')
end

list = {'Separation','Period'};
[indx,tf] = listdlg('PromptString','Choose what to plot in the ordinate',...
                    'SelectionMode','single',...
                    'ListString',list);

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

% Evaluate functions
separation_AU = linspace(0.01,11,num_dots);
m3_Mardling = logspace(0.0,3,num_dots);
temp_Mardling = aout_crit_Mardling(separation_inner_AU,Minner,m3_Mardling);

if indx==1
    y_Mardling = temp_Mardling;
elseif indx==2
    y_Mardling = orbital_period_yr(temp_Mardling,Minner+m3_Mardling).*AstroConstants.yr_to_d;
else
    display('Faulty choice.')
end

% Find systems of interest (so that the plot is not crowded over the edges)
indexOfInterest = find((m3>1 & m3<10^4) & (aout<=10));
% limit_indices = find((fraction>0.1)&(fraction<0.2)&(m3>4.0));
eta_accuracy = 0.02;
ioi_ETA_unity = find(ETA<1+eta_accuracy & ETA>1-eta_accuracy);

periastron_accuracy_Rsun = 0.1;
ioi_contact = find(0.5002.*periastron_min_Rsun<=radius_Rsun+periastron_accuracy_Rsun & 0.5002.*periastron_min_Rsun>=radius_Rsun-periastron_accuracy_Rsun);

% Plot
solar=char(9737);
sz = 140.0;
lw=3.0;
fs=18;

lines1 = [         0    0.4470    0.7410];
lines2 = [    0.8500    0.3250    0.0980];
lines3 = [    0.9290    0.6940    0.1250];
lines4 = [    0.4940    0.1840    0.5560];
lines5 = [    0.4660    0.6740    0.1880];
lines6 = [    0.3010    0.7450    0.9330];
lines7 = [    0.6350    0.0780    0.1840];
grey   = 0.5.*[1 1 1];

pos_hor_6_d = 9.5;
pos_ver_6_d = 0.41;
pos_hor_30_d = 3.5;
pos_ver_30_d = 1.1;

pos_hor     = 70;
pos_wide    = 6.0;
pos_close   = 0.5;
pos_unst    = 0.15;

xLims       = [1 220];
xLimsLabels = {'1','10','100'};

if indx==1
    yLims       = [0.09 11];
    % yLimsLabels = {'0.1','1','10'};
elseif indx==2
    yLims       = [2 1000];
    % yLimsLabels = {'0.1','1','10'};
else
    display('Plotting a faulty choice.')
end


clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
% t = tiledlayout(2,3,'TileSpacing','Compact');
t = tiledlayout(2,3);

% 1. e_max
nexttile
hold on
box on
cbar1 = colorbar;
cbar1.FontSize = fs;
cbar1.TickLabelInterpreter = 'latex';
cbar1.Location = "northoutside";
cbar1.Label.Interpreter = 'latex';
cbar1.Label.String = '$e_{max}$';
colormap(flip(pink))

ax1 = gca;
ax1.FontSize = fs;
ax1.XLim = xLims;
ax1.YLim = yLims;
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.XTickLabels = xLimsLabels;
% ax1.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,e_max(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled')
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,e_max(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled')
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)
legend( '$\eta \approx 1$',...
        '$M_{1}=M_{2}=M_{3}$',...
        '$a_{\rm{crit}}$',...
        'interpreter','latex',...
        'location','West')

% 2. ZKL timescale
nexttile
hold on
box on
cbar2 = colorbar;
cbar2.FontSize = fs;
cbar2.TickLabelInterpreter = 'latex';
cbar2.Location = "northoutside";
cbar2.Label.Interpreter = 'latex';
cbar2.Label.String = '$\log_{10} (\tau_{\rm{ZLK}}/\rm{yr})$';

colormap(flip(pink))

ax2 = gca;
ax2.FontSize = fs;
ax2.XLim = xLims;
ax2.YLim = yLims;
ax2.XScale = 'log';
ax2.YScale = 'log';
ax2.XTickLabels = xLimsLabels;
% ax2.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,log10(tau_ZKL(indexOfInterest)),'s','Filled')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,log10(tau_ZKL(indexOfInterest)),'s','Filled')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)

% 3. eta
nexttile
hold on
box on
cbar3 = colorbar;
cbar3.FontSize = fs;
cbar3.TickLabelInterpreter = 'latex';
cbar3.Location = "northoutside";
cbar3.Label.Interpreter = 'latex';
cbar3.Label.String = '$\log_{10} \eta$';

colormap(flip(pink))

ax3 = gca;
ax3.FontSize = fs;
ax3.XLim = xLims;
ax3.YLim = yLims;
ax3.XScale = 'log';
ax3.YScale = 'log';
ax3.XTickLabels = xLimsLabels;
% ax3.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,log10(ETA(indexOfInterest)),'s','Filled')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,log10(ETA(indexOfInterest)),'s','Filled')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)

% 4. f_contact
nexttile
hold on
box on
cbar4 = colorbar;
cbar4.FontSize = fs;
cbar4.TickLabelInterpreter = 'latex';
cbar4.Location = "northoutside";
cbar4.Label.Interpreter = 'latex';
cbar4.Label.String = '$f_{\rm{contact}}$';

colormap(flip(pink))

ax4 = gca;
ax4.FontSize = fs;
ax4.XLim = xLims;
ax4.YLim = yLims;
ax4.XScale = 'log';
ax4.YScale = 'log';
ax4.XTickLabels = xLimsLabels;
% ax4.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,fraction(indexOfInterest),'s','Filled')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled')
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,fraction(indexOfInterest),'s','Filled')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)

% 5. rps
nexttile
hold on
box on
cbar5 = colorbar;
cbar5.FontSize = fs;
cbar5.TickLabelInterpreter = 'latex';
cbar5.Location = "northoutside";
cbar5.Label.Interpreter = 'latex';
cbar5.Label.String = '$r_{\rm{p,min}}=a_{\rm{in}}(1-e_{\rm{max}})$';

colormap(flip(pink))

ax5 = gca;
ax5.FontSize = fs;
ax5.XLim = xLims;
ax5.YLim = yLims;
ax5.XScale = 'log';
ax5.YScale = 'log';
ax5.XTickLabels = xLimsLabels;
% ax5.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,periastron_min_Rsun(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled','HandleVisibility','off')    
    scatter(m3(ioi_contact),aout(ioi_contact),sz,'g','s','Filled')    
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,periastron_min_Rsun(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled','HandleVisibility','off')    
    scatter(m3(ioi_contact),orbital_period_outer_days(ioi_contact),sz,'g','s','Filled')      
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)
legend( '$R_{*} \approx r_{\rm{p}}/2$',...
        'interpreter','latex',...
        'location','West')

% 6. cos_inc
nexttile
hold on
box on
cbar6 = colorbar;
cbar6.FontSize = fs;
cbar6.TickLabelInterpreter = 'latex';
cbar6.Location = "northoutside";
cbar6.Label.Interpreter = 'latex';
cbar6.Label.String = '$\cos(i)|_{e_{\rm{max}}}$';
colormap(flip(pink))

ax6 = gca;
ax6.FontSize = fs;
ax6.XLim = xLims;
ax6.YLim = yLims;
ax6.XScale = 'log';
ax6.YScale = 'log';
ax6.XTickLabels = xLimsLabels;
% ax6.YTickLabels = yLimsLabels;

scatter(1000,1000,1,0,'HandleVisibility','off')
scatter(1000,1000,1,1,'HandleVisibility','off')
xlabel('$M_3/M_\odot$','Interpreter','latex','FontSize',fs)

if indx==1
    scatter(m3(indexOfInterest),aout(indexOfInterest),sz,cos_inc(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),aout(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$a_{\rm{out}}/\rm{AU}$','Interpreter','latex','FontSize',fs)
elseif indx==2
    scatter(m3(indexOfInterest),orbital_period_outer_days(indexOfInterest),sz,cos_inc(indexOfInterest),'s','Filled','HandleVisibility','off')
    scatter(m3(ioi_ETA_unity),orbital_period_outer_days(ioi_ETA_unity),sz,grey,'s','Filled')    
    ylabel('$P_{\rm{out}}/\rm{d}$','Interpreter','latex','FontSize',fs)    
else
    display('Plotting a faulty choice.')
end

xline(mass_Msun,'r--','LineWidth',lw)
plot(m3_Mardling, y_Mardling,'-','LineWidth',lw,'Color',lines6)
legend( '$R_{*} \approx r_{\rm{p}}/2$',...
        'interpreter','latex',...
        'location','West')

% Save
if save_flag
    print(gcf,plot_label,'-dpng','-r250');
end

end
