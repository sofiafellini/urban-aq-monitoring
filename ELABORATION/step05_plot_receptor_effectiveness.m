% =========================================================================
% PLOT: INDIVIDUAL NODE PERFORMANCE MAP
% =========================================================================
%
% This script visualises the individual node performance score R_a on the
% urban road network graph. Each node is coloured according to its R_a
% value, which summarises how useful that node is as a receptor across all
% pairs or triplets.
%
% A higher R_a value (warmer colour) indicates that placing a sensor at
% that node tends to yield better source identification performance.
%
% -------------------------------------------------------------------------
% INPUTS:
%   ../INPUT/<...>.mat                              - city graphs
%   single_node_<city>_<dir>_delta0<d>.txt          - individual node scores
%                                                     (output of step01_receptor_pairs.m)
%   triplets_single_node_<city>_<dir>_delta0<d>.txt - individual node scores
%                                                     (output of step03_receptor_triplets.m)
%
% =========================================================================

clear; clc; close all;

addpath('..')
addpath('PLOT_FUNCTIONS')

% -------------------------------------------------------------------------
% USER PARAMETERS
% -------------------------------------------------------------------------

city_names   = {'Firenze', 'Parigi', 'New_York', 'Lione'};
phi_rotation = [34, 23, 124, 169];    % rotation angle of main street axis [deg]
wind_dirs    = [360, 45, 90, 135, 180, 225, 270, 315];  % wind directions [deg]
delta        = 0.1;

% Select which R_a scores to visualise:
%   'pairs'    - individual node performance summed over all receptor PAIRS
%                (output of step01_receptor_pairs.m)
%   'triplets' - individual node performance summed over all receptor TRIPLETS
%                (output of step03_receptor_triplets.m)
receptor_mode = 'pairs';

% Colour map
my_map_hsv = customcolormap_preset('pasteljet');

% View angles for each city [azimuth, elevation]
rot_angles = [-15, -20, -10,  0];
elev_angles = [ 80,  85,  90, 90];

% -------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------

for idx_city = 1 %choose the city in city_names

    for idx_dir = 7 %choose the wind direction in wind_dirs

        f1 = figure('PaperPosition', [0 0 50 20]);

        % Effective wind direction
        wind_dir_effective = wind_dirs(idx_dir) + phi_rotation(idx_city);
        if wind_dir_effective > 360
            wind_dir_effective = wind_dir_effective - 360;
        end

        % -----------------------------------------------------------------
        % LOAD INDIVIDUAL NODE PERFORMANCE SCORES
        % -----------------------------------------------------------------
        if strcmp(receptor_mode, 'triplets')
            file_nodes = ['triplets_single_node_' city_names{idx_city} '_' ...
                num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
        else
            file_nodes = ['single_node_' city_names{idx_city} '_' ...
                num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
        end

        node_data = table2array(readtable(file_nodes));
        node_scores = node_data(:, 2);   % R_a values

        % -----------------------------------------------------------------
        % LOAD URBAN DISPERSION DATA (for graph structure and node coords)
        % -----------------------------------------------------------------
        input_file = fullfile('..', 'INPUT', ...
            ['5_', num2str(wind_dir_effective), '_C0_Cth_10_', ...
             city_names{idx_city}, '_25_03.mat']);
        load(input_file)

        % -----------------------------------------------------------------
        % PLOT GRAPH WITH NODE COLOUR = R_a
        % -----------------------------------------------------------------
        figure(f1)
        h = plot(G_str, ...
            'XData',      X_node_a, ...
            'YData',      Y_node_a, ...
            'NodeColor',  'k', ...
            'MarkerSize', 3.5, ...
            'ArrowSize',  0.01);

        % Update the figure title
        title(gca,['Performance of nodes as receptors: ' city_names{idx_city} ', $\Phi=$ ' num2str(wind_dirs(idx_dir)) '$^\circ$'],'FontSize', 12);

        % Colour each node by its performance score R_a
        h.NodeCData = node_scores;
        colormap(my_map_hsv)
        view(rot_angles(idx_city), elev_angles(idx_city))
        % Remove axes ticks and labels
        set(gca, 'XTickLabel', '', 'XTick', '');
        set(gca, 'YTickLabel', '', 'YTick', '');
        %set(gca,'visible','off')

        % Colorbar
        figure(f1)
        cb = colorbar('southoutside');
        title(cb, '$\mathcal{R}_a$', 'Interpreter', 'Latex', 'FontSize', 20)
    end
end