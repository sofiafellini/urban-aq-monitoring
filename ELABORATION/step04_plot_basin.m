% =========================================================================
% PLOT: DETECTION BASIN OF THE BEST RECEPTOR PAIR (OR TRIPLET)
% =========================================================================
%
% This script visualises the urban road network graph and highlights the
% detection basin of the best receptor pair identified by step01_receptor_pairs.m.
%
% It produces the figures showing the optimal receptor pair and its
% detection basin for each city, as presented in Figures 3 and 4 of:
%
%   Carlini, B., Salizzoni, P., Ridolfi, L., Fellini, S. (2026).
%   A model-independent matrix method for optimizing urban air quality
%   monitoring networks. Environmental Modelling & Software.
%
% For the best recetpor pair, the code highlights nodes as follows (as in Figure 4 of the paper):
%   - Dark grey  : nodes outside the detection basin
%   - Blue       : nodes with Pi = 0 (no source identification capability)
%   - Yellow     : nodes with Pi > 0 (partial source identification)
%   - Green      : nodes with Pi = 1 (full, unambiguous source identification)
%   - Dark red   : receptor nodes (sensor locations)
%
% For the best recetpor triplet, the code highlights nodes as follows (Figure 8 of the paper):
% the detection basin is split into three sub-basins (ab, bc, ac), each shown in a different colour.
%
% -------------------------------------------------------------------------
% NOTE ON DIRECTORY STRUCTURE
%
% This script must be run from the ELABORATION\ folder.
% Input .mat files are loaded from ..\INPUT\
%
% -------------------------------------------------------------------------
% INPUTS:
%   ..\INPUT\5_<dir>_C0_Cth_10_<city>_25_03.mat  - urban network data
%                                                   (graph structure and
%                                                   node coordinates)
%   M_bestpair_<city>_<dir>_delta0<d>.txt         - metric matrix sorted
%                                                   by R_ab (output of
%                                                   step01_receptor_pairs.m)
%   v_certain_nodes_<city>_<dir>_delta0<d>.txt    - certain source nodes
%                                                   (output of
%                                                   step01_receptor_pairs.m)
%   v_uncertain_<city>_<dir>_delta0<d>.txt        - uncertain source nodes
%                                                   (output of
%                                                   step01_receptor_pairs.m)
%   v_zero_nodes_<city>_<dir>_delta0<d>.txt       - basin nodes with Pi = 0
%                                                   (output of
%                                                   step01_receptor_pairs.m)
%
%   (for triplet visualisation)
%   M_triplets_<city>_<dir>_delta0<d>.txt         - metric matrix for
%                                                   triplets (output of
%                                                   step03_receptor_triplets.m)
%   best_triplet_basin_<city>_<dir>_delta0<d>.txt - detection basins of
%                                                   the best triplet
%                                                   (output of
%                                                   step03_receptor_triplets.m)
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
delta        = 0.1;                   % relative measurement error [-]

% Select which basin to visualise:
%   'pairs'    - basin of the best receptor pair
%                (output of step01_receptor_pairs.m)
%   'triplets' - basin of the best receptor triplet
%                (output of step03_receptor_triplets.m)
receptor_mode = 'pairs';

% -------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------

for idx_city = 1 %choose the city in city_names

    for idx_dir = 7 %choose the wind direction in wind_dirs

        % Effective wind direction (accounting for street-axis rotation)
        wind_dir_effective = wind_dirs(idx_dir) + phi_rotation(idx_city);
        if wind_dir_effective > 360
            wind_dir_effective = wind_dir_effective - 360;
        end

        % -----------------------------------------------------------------
        % LOAD URBAN NETWORK DATA (graph structure and node coordinates)
        % -----------------------------------------------------------------
        input_file = fullfile('..', 'INPUT', ...
            ['5_', num2str(wind_dir_effective), '_C0_Cth_10_', ...
             city_names{idx_city}, '_25_03.mat']);
        load(input_file)

        % -----------------------------------------------------------------
        % PLOT BASE GRAPH
        % -----------------------------------------------------------------
        figure

        h = plot(G_str, ...
            'XData',       X_node_a, ...
            'YData',       Y_node_a, ...
            'NodeColor',   [0.3 0.3 0.3], ...
            'MarkerSize',  2, ...
            'ShowArrows',  'off');

        if strcmp(receptor_mode, 'pairs')

            % -----------------------------------------------------------------
            % LOAD PAIR RESULTS AND HIGHLIGHT DETECTION BASIN
            % -----------------------------------------------------------------
    
            % Build dynamic filenames consistent with step01_receptor_pairs.m outputs
            file_bestpair   = ['M_bestpair_'      city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
            file_certain    = ['v_certain_nodes_' city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
            file_uncertain  = ['v_uncertain_'     city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
            file_zero       = ['v_zero_nodes_'    city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
    
            M_bestpair      = table2array(readtable(file_bestpair));
            certain_nodes   = table2array(readtable(file_certain));
            uncertain_nodes = table2array(readtable(file_uncertain));
            zero_nodes      = table2array(readtable(file_zero));
    
            % Keep only valid (positive) node indices
            certain_nodes   = certain_nodes(certain_nodes > 0);
            uncertain_nodes = uncertain_nodes(uncertain_nodes > 0);
            zero_nodes      = zero_nodes(zero_nodes > 0);
    
            % Colour nodes by source identification capability (Figures 3-4):
            %   blue     = Pi = 0  (in basin, not discriminable)
            %   yellow   = Pi > 0  (uncertain nodes)
            %   green    = Pi = 1  (certain nodes)
            %   dark red = receptor locations
            highlight(h, zero_nodes,      'NodeColor', [0.3010 0.7450 0.9330], 'MarkerSize', 3.2)
            highlight(h, uncertain_nodes, 'NodeColor', [0.9290 0.6940 0.1250], 'MarkerSize', 3.2)
            highlight(h, certain_nodes,   'NodeColor', [0.4660 0.6740 0.1880], 'MarkerSize', 3.2)
            highlight(h, [M_bestpair(end, 1), M_bestpair(end, 2)], ...
                'NodeColor', [0.6350 0.0780 0.1840], 'MarkerSize', 4)
    
    
            % --- Legend Construction Block (Scientific Publication Standard) ---
            hold on;
    
            % Create invisible dummy points for the legend using your specific RGB colors
            h_zero   = scatter(NaN, NaN, 30, [0.3010 0.7450 0.9330], 'filled'); % Light blue
            h_uncert = scatter(NaN, NaN, 30, [0.9290 0.6940 0.1250], 'filled'); % Gold/Yellow
            h_cert   = scatter(NaN, NaN, 30, [0.4660 0.6740 0.1880], 'filled'); % Green
            h_recept = scatter(NaN, NaN, 50, [0.6350 0.0780 0.1840], 'filled', 'Marker', 's'); % Maroon square
    
            % Legend labels mapped to your manuscript's mathematical definitions
            legend_entries = { ...
                'Nodes in detection basin ($P_i=0$)', ...
                'Uncertain source locations ($P_i>0$)', ...
                'Reliable source locations (($P_i=1$))', ...
                'Optimal receptor pair'};
    
            % Generate the legend with professional formatting
            lgd = legend([h_zero, h_uncert, h_cert, h_recept], legend_entries, ...
                'Interpreter','Latex','Location', 'bestoutside', 'FontSize', 10);
            title(lgd, 'Source Identification Mapping');

        else 

            % -----------------------------------------------------------------
            % TRIPLET VISUALISATION (uncomment to produce Figure 8)
            % Requires output files from step03_receptor_triplets.m
            % -----------------------------------------------------------------
            file_triplets = ['M_triplets_'         city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
            file_basins   = ['best_triplet_basin_' city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta*10) '.txt'];
    
            M_triplets = table2array(readtable(file_triplets));
            v_basins   = table2array(readtable(file_basins));
            v_ab = v_basins(:, 1);
            v_bc = v_basins(:, 2);
            v_ac = v_basins(:, 3);
    
            % Sub-basin ab (light blue), bc (green), ac (yellow)
            highlight(h, v_ab(~isnan(v_ab)), 'NodeColor', [0.4  0.7  1.0 ], 'MarkerSize', 3.2)
            highlight(h, v_bc(~isnan(v_bc)), 'NodeColor', [0.4660 0.6740 0.1880], 'MarkerSize', 3.2)
            highlight(h, v_ac(~isnan(v_ac)), 'NodeColor', [0.9290 0.6940 0.1250], 'MarkerSize', 3.2)
            highlight(h, [M_triplets(end,1), M_triplets(end,2), M_triplets(end,3)], ...
                'NodeColor', [0.6350 0.0780 0.1840], 'MarkerSize', 4)

             % --- Legend Construction Block (Scientific Publication Standard) ---
            hold on;
    
            % Create invisible dummy points for the legend using your specific RGB colors
            h_ab   = scatter(NaN, NaN, 30, [0.4  0.7  1.0 ], 'filled'); % Light blue
            h_bc = scatter(NaN, NaN, 30, [0.4660 0.6740 0.1880], 'filled'); % Gold/Yellow
            h_ac   = scatter(NaN, NaN, 30, [0.9290 0.6940 0.1250], 'filled'); % Green
            h_recept = scatter(NaN, NaN, 50, [0.6350 0.0780 0.1840], 'filled', 'Marker', 's'); % Maroon square
    
            % Legend labels mapped to your manuscript's mathematical definitions
            legend_entries = { ...
                'Basin ab', ...
                'Basin bc', ...
                'Basin ac', ...
                'Optimal receptor triplet'};
    
            % Generate the legend with professional formatting
            lgd = legend([h_ab, h_bc, h_ac, h_recept], legend_entries, ...
                'Interpreter','Latex','Location', 'bestoutside', 'FontSize', 10);
            title(lgd, 'Source Identification Mapping');

       
        end

       % Update the figure title
       title(gca,['Detection Basin Analysis: ' city_names{idx_city} ', $\Phi=$ ' num2str(wind_dirs(idx_dir)) '$^\circ$'],'FontSize', 12);
        % --- End Legend Block ---
    
    
        % Remove axes ticks and labels
        set(gca, 'XTickLabel', '', 'XTick', '');
        set(gca, 'YTickLabel', '', 'YTick', '');
    
    end

end