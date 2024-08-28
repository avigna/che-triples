function plot_initial_distributions(debug_flag)
tic;

% DATA
% Load adata (provided by Abinaya)
% Binaries
filename_binary     = '../data/initial_distributions/inipop_binary.hdf5';
bin_AP1s            = h5read(filename_binary,'/AP1s');
bin_INCL1s          = h5read(filename_binary,'/INCL1s');
bin_LAN1s           = h5read(filename_binary,'/LAN1s');
bin_a1s             = h5read(filename_binary,'/a1s');
bin_e1s             = h5read(filename_binary,'/e1s');
bin_m1s             = h5read(filename_binary,'/m1s');
bin_m2s             = h5read(filename_binary,'/m2s');

% Triples
filename_triple     = '../data/initial_distributions/inipop_triple.hdf5';
tri_AP1s            = h5read(filename_triple,'/AP1s');
tri_AP2s            = h5read(filename_triple,'/AP2s');
tri_INCL1s          = h5read(filename_triple,'/INCL1s');
tri_INCL2s          = h5read(filename_triple,'/INCL2s');
tri_LAN1s           = h5read(filename_triple,'/LAN1s');
tri_LAN2s           = h5read(filename_triple,'/LAN2s');
tri_a1s             = h5read(filename_triple,'/a1s');
tri_a2s             = h5read(filename_triple,'/a2s');
tri_e1s             = h5read(filename_triple,'/e1s');
tri_e2s             = h5read(filename_triple,'/e2s');
tri_m1s             = h5read(filename_triple,'/m1s');
tri_m2s             = h5read(filename_triple,'/m2s');
tri_m3s             = h5read(filename_triple,'/m3s');

% Quadruples
filename_quadruple  = '../data/initial_distributions/inipop_quadruple.hdf5';
quad_AP1s           = h5read(filename_quadruple,'/AP1s');
quad_AP2s           = h5read(filename_quadruple,'/AP2s');
quad_AP3s           = h5read(filename_quadruple,'/AP3s');
quad_INCL1s         = h5read(filename_quadruple,'/INCL1s');
quad_INCL2s         = h5read(filename_quadruple,'/INCL2s');
quad_INCL3s         = h5read(filename_quadruple,'/INCL3s');
quad_LAN1s          = h5read(filename_quadruple,'/LAN1s');
quad_LAN2s          = h5read(filename_quadruple,'/LAN2s');
quad_LAN3s          = h5read(filename_quadruple,'/LAN3s');
quad_a1s            = h5read(filename_quadruple,'/a1s');
quad_a2s            = h5read(filename_quadruple,'/a2s');
quad_a3s            = h5read(filename_quadruple,'/a3s');
quad_e1s            = h5read(filename_quadruple,'/e1s');
quad_e2s            = h5read(filename_quadruple,'/e2s');
quad_e3s            = h5read(filename_quadruple,'/e3s');
quad_m1s            = h5read(filename_quadruple,'/m1s');
quad_m2s            = h5read(filename_quadruple,'/m2s');
quad_m3s            = h5read(filename_quadruple,'/m3s');
quad_m4s            = h5read(filename_quadruple,'/m4s');

% Subselections
ioi_compact_trips = find(tri_a1s<=1 & tri_a2s<=100);
ioi_compact_quads = find(quad_a1s<=1 & quad_a3s<=100);
% nnz(find(tri_a1s<1))
% 
% nnz(find(quad_a1s<1))

% Count systems
if debug_flag
    length(bin_AP1s)
    length(tri_AP1s)
    length(quad_AP1s)
    length(ioi_compact_trips)
    length(ioi_compact_quads)    
end

% Create CDFs
[a_in_bin, a_in_bin_CDF]    = create_empirical_CDF(bin_a1s, ones(size(bin_a1s)));
[a_in_tri, a_in_tri_CDF]    = create_empirical_CDF(tri_a1s, ones(size(tri_a1s)));
[a_in_quad, a_in_quad_CDF]  = create_empirical_CDF(quad_a1s, ones(size(quad_a1s)));

[a_out_tri, a_out_tri_CDF]  = create_empirical_CDF(tri_a2s, ones(size(tri_a2s)));
[a_out_quad, a_out_quad_CDF]= create_empirical_CDF(quad_a3s, ones(size(quad_a3s)));

% PLOT
lw=2.0;
fs=14;
sz=2.0;

lines1 = [         0    0.4470    0.7410];
lines2 = [    0.8500    0.3250    0.0980];
lines3 = [    0.9290    0.6940    0.1250];
lines4 = [    0.4940    0.1840    0.5560];
lines5 = [    0.4660    0.6740    0.1880];
lines6 = [    0.3010    0.7450    0.9330];
lines7 = [    0.6350    0.0780    0.1840];

clf
set(gca,'FontName','Helvetica')
box on
hold on
plot(log10(a_in_bin),a_in_bin_CDF,'LineWidth',lw,'Color',lines1)
plot(log10(a_in_tri),a_in_tri_CDF,'LineWidth',lw,'Color',lines2)
plot(log10(a_in_quad),a_in_quad_CDF,'LineWidth',lw,'Color',lines3)

plot(log10(a_out_tri),a_out_tri_CDF,'--','LineWidth',lw,'Color',lines2)
plot(log10(a_out_quad),a_out_quad_CDF,'--','LineWidth',lw,'Color',lines3)

xlabel('log_{10} a/au','FontSize',fs)
ylabel('CDF','FontSize',fs)

legend( 'a_{bin}',...
        'a_{tri,in}',...
        'a_{quad,in}',...
        'a_{tri,out}',...
        'a_{quad,out}',...
        'Location','SouthEast',...
        'FontSize',fs)

ax1 = gca;
ax1.FontName = 'Helvetica';
ax1.FontSize = fs;
ax1.XLim = [-1 4.5];
ax1.XTick = [-1:4];
ax1.XGrid = 'on';
ax1.YGrid = 'off';
ax1.XAxis.FontName= 'Helvetica';
ax1.YAxis.FontName= 'Helvetica';

figure()
clf

scatter(log10(tri_a1s(ioi_compact_trips)),log10(tri_a2s(ioi_compact_trips)),sz,log10(tri_m1s(ioi_compact_trips)))
scatter(log10(quad_a1s(ioi_compact_quads)),log10(quad_a2s(ioi_compact_quads)),sz,log10(quad_m1s(ioi_compact_quads)))


colorbar



toc;
end