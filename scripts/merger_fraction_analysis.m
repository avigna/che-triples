function merger_fraction_analysis(mass, separation_inner, debug_flag, save_flag)
% Function created by Alejandro Vigna-Gomez
% [mass] = Msun
% [separation_inner] = Rsun
% debug_flag: bool
% save_flag: bool

% From COMPAS (and used in Bin's simulation) a55_au = 20.16/auToRsun
% merger_fraction_analysis(55,20.16,0,0)

% Macros
year_to_day = 365.24;
auToRsun = 215.032;

% Functions
aout_crit_Mardling = @(ain,Minner,M3)  2.8.*ain.*(1+(M3./Minner)).^(2.0/5);

% From Kepler's law
isoMass = @(a_out,P_out,M_bin) ((a_out.^3)*(P_out.^-2))-M_bin;

% calculateRocheRadius(Md,Ma)
%   Roche radius following Eggleton 1983.
%   Roche radius as seen by star 1, e.g. donor star.
%   q = Md/Ma
% Calculate Roche radius
% roche_radius_inner_binary   = @(a_out,M_bin,M_out) a_out.*calculateRocheRadius(M_bin,M_out);
% roche_radius_outer_star     = @(a_out,M_bin,M_out) a_out.*calculateRocheRadius(M_out,M_bin);

% Load data
filename=strcat('../data/Dynamics/m1_',num2str(mass),'.txt');
M=load(filename);
m3 = M(:,2);
aout = M(:,3);
tauKL = M(:,4);
fraction =  M(:,5);

% Print values
if debug_flag
    filename
    fprintf('Mass_1 = Mass 2 = %f',mass)
    fprintf('\n')
    fprintf('a_{inner} = %f', separation_inner)
    fprintf('\n')
    fprintf('min(f_merge) = %f',min(fraction))
    fprintf('\n')
    fprintf('max(f_merge) = %f',max(fraction))
    fprintf('\n')
end

% Some quantities
num_dots = 1000;
separation_inner_AU = separation_inner/auToRsun;
Minner = 2*mass;
Pdays = [1 5 6 10 30 50 100 500 1000];
Pyears = Pdays./year_to_day;

if debug_flag
    Pdays(2)
    Pdays(3)
    Pdays(5)
end

% Evaluate functions
separation_AU = linspace(0.01,11,num_dots);
m3_Mardling = logspace(0.0,3,num_dots);
y_Mardling = aout_crit_Mardling(separation_inner_AU,Minner,m3_Mardling);


% Find systems of interest (so that the plot is not crowded over the edges)
indexOfInterest = find((m3>1 & m3<10^4) & (aout<=10));

limit_indices = find((fraction>0.1)&(fraction<0.2)&(m3>4.0));
f = fit(m3(limit_indices),aout(limit_indices),'b*x^m');

aout(limit_indices)
fit_separation_AU   = f.b.*m3(limit_indices).^(f.m)
fit_m3              = m3(limit_indices);
fit_period_yr       = sqrt((aout(limit_indices).^3)./(m3(limit_indices)+Minner));
fit_period_d        = fit_period_yr.*year_to_day

round(m3(limit_indices))
round(fit_period_d)
% 
% f2 = fit(m3(limit_indices),fit_period_d,'b*x^m')
% 
% fit_period_d
% f2.b.*m3(limit_indices).^(f2.m)

% Find intersection stability and tidal radius
% valueFromPlot = 0.1;
% indexOfIntersection = find(x_tidal_radius>valueFromPlot);
% Bin_aout_crit = aout_crit_Mardling(separation_inner_AU,Minner,m3);
% fraction(find(aout<Bin_aout_crit-0.01))=0.0;

% Plot
solar=char(9737);
sz = 80.0;
lw=2.0;
fs=14;

lines1 = [         0    0.4470    0.7410];
lines2 = [    0.8500    0.3250    0.0980];
lines3 = [    0.9290    0.6940    0.1250];
lines4 = [    0.4940    0.1840    0.5560];
lines5 = [    0.4660    0.6740    0.1880];
lines6 = [    0.3010    0.7450    0.9330];
lines7 = [    0.6350    0.0780    0.1840];

pos_hor_6_d = 9.5;
pos_ver_6_d = 0.41;
pos_hor_30_d = 3.0;
pos_ver_30_d = 1.2;

pos_hor     = 70;
pos_wide    = 6.0;
pos_close   = 0.65;
pos_unst    = 0.17;

clf
t = tiledlayout(3,1,'TileSpacing','Compact');
nexttile
hold on
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink))

y_Mardling = [0.1 y_Mardling];
m3_Mardling = [1 m3_Mardling];
unstableRegion_Mardling = patch([m3_Mardling 10^4*ones(size(m3_Mardling))],[y_Mardling fliplr(y_Mardling)],'k');
unstableRegion_Mardling.EdgeColor = 'none';
unstableRegion_Mardling.FaceColor = 0.5.*[1 1 1];

% plot(isoMass(separation_AU, Pyears(2), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% plot(isoMass(separation_AU, Pyears(4), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
plot(isoMass(separation_AU, Pyears(5), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% text(pos_hor_6_d,pos_ver_6_d,'P_{out} = 6 d','Color',lines6,'Fontsize',fs)
text(pos_hor_30_d,pos_ver_30_d,'P_{out} = 30 d','Color',lines6,'Fontsize',fs)

% xlabel(['M_3/M_',solar])
ylabel('a_{out}/AU')
title(['M_1=M_2=',num2str(mass),' M_',solar,' at Z=0.0001']);

text(pos_hor,pos_wide,['Wide',newline,'Triples'],'Color','k','Fontsize',fs)
text(pos_hor,pos_close,['Close',newline,'Triples'],'Color','k','Fontsize',fs)
text(pos_hor,pos_unst,['Dynamically',newline,'Unstable'],'Color','w','Fontsize',fs)

ax1 = gca;
ax1.FontSize = fs;
ax1.XLim = [1 200];
ax1.YLim = [0.1 10];
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.XTickLabels = {'1','10','100'};
ax1.YTickLabels = {'0.1','1','10'};
box on

nexttile
hold on
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink))

scatter(m3(indexOfInterest),aout(indexOfInterest),sz,fraction(indexOfInterest),'s','Filled')
cbar.Label.String = 'f_{contact}';

y_Mardling = [0.1 y_Mardling];
m3_Mardling = [1 m3_Mardling];
unstableRegion_Mardling = patch([m3_Mardling 10^4*ones(size(m3_Mardling))],[y_Mardling fliplr(y_Mardling)],'k');
unstableRegion_Mardling.EdgeColor = 'none';
unstableRegion_Mardling.FaceColor = 0.5.*[1 1 1];

% plot(isoMass(separation_AU, Pyears(2), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% plot(isoMass(separation_AU, Pyears(4), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
plot(isoMass(separation_AU, Pyears(5), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% text(pos_hor_6_d,pos_ver_6_d,'P_{out} = 6 d','Color',lines6,'Fontsize',fs)
text(pos_hor_30_d,pos_ver_30_d,'P_{out} = 30 d','Color',lines6,'Fontsize',fs)


% plot(m3(limit_indices),aout(limit_indices),'k')
% xlabel(['M_3/M_',solar])
ylabel('a_{out}/AU')

text(pos_hor,pos_wide,['Wide',newline,'Triples'],'Color','k','Fontsize',fs)
text(pos_hor,pos_close,['Close',newline,'Triples'],'Color','w','Fontsize',fs)
text(pos_hor,pos_unst,['Dynamically',newline,'Unstable'],'Color','w','Fontsize',fs)

% if debug_flag
%     plot(m3(limit_indices),f.b.*m3(limit_indices).^(f.m),'k')
% end

ax1 = gca;
ax1.FontSize = fs;
ax1.XLim = [1 200];
ax1.YLim = [0.1 10];
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.XTickLabels = {'1','10','100'};
ax1.YTickLabels = {'0.1','1','10'};
box on

nexttile
hold on
cbar = colorbar;
cbar.FontSize = fs;
colormap(flip(pink))

scatter(m3(indexOfInterest),aout(indexOfInterest),sz,log10(tauKL(indexOfInterest)),'s','Filled')
cbar.Label.String = 'log_{10} (\tau_{ZLK}/yr)';

y_Mardling = [0.1 y_Mardling];
m3_Mardling = [1 m3_Mardling];
unstableRegion_Mardling = patch([m3_Mardling 10^4*ones(size(m3_Mardling))],[y_Mardling fliplr(y_Mardling)],'k');
unstableRegion_Mardling.EdgeColor = 'none';
unstableRegion_Mardling.FaceColor = 0.5.*[1 1 1];

% plot(isoMass(separation_AU, Pyears(2), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% plot(isoMass(separation_AU, Pyears(4), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
plot(isoMass(separation_AU, Pyears(5), Minner),separation_AU,'-','LineWidth',lw,'Color',lines6)
% text(pos_hor_6_d,pos_ver_6_d,'P_{out} = 6 d','Color',lines6,'Fontsize',fs)
text(pos_hor_30_d,pos_ver_30_d,'P_{out} = 30 d','Color',lines6,'Fontsize',fs)

xlabel(['M_3/M_',solar])
ylabel('a_{out}/AU')

text(pos_hor,pos_wide,['Wide',newline,'Triples'],'Color','k','Fontsize',fs)
text(pos_hor,pos_close,['Close',newline,'Triples'],'Color','k','Fontsize',fs)
text(pos_hor,pos_unst,['Dynamically',newline,'Unstable'],'Color','w','Fontsize',fs)

ax1 = gca;
ax1.FontSize = fs;
ax1.XLim = [1 200];
ax1.YLim = [0.1 10];
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.XTickLabels = {'1','10','100'};
ax1.YTickLabels = {'0.1','1','10'};
box on

% Save
plot_label = '../plots/test.png';
if save_flag
    print(gcf,plot_label,'-dpng','-r250');
end

end