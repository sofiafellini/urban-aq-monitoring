% =========================================================================
% RECEPTOR TRIPLET IDENTIFICATION FOR URBAN POLLUTANT SOURCE DETECTION
% =========================================================================
%
% This script is the fourth and final step in the receptor identification
% pipeline, and must be run after:
%   1. step01_receptor_pairs.m
%   2. step02_receptor_pairs_basin.m
%
% For each candidate receptor triplet (a, b, c), the combined detection
% basin is assembled by merging the basins of the three constituent pairs
% (ab, ac, bc), resolving overlaps in favour of the pair with the highest
% combined performance index R.
%
% The combined performance index of the triplet is defined as R_abc
%   
%
% The best triplet (highest R_abc) is then characterised in detail, and
% its detection basin is partitioned into three non-overlapping sub-basins,
% one per constituent pair.
%
% -------------------------------------------------------------------------
% INPUTS:
%   M_<city>_<dir>_delta<delta*10>_original.txt    - pair metric matrix (unsorted),
%                                            output of step01_receptor_pairs.m
%                                            cols: [a, b, B_(ab,n), P_ab, R_ab]
%   pairs_<city>_<dir>_basin_nodes.txt     - basin node lists per pair,
%                                            output of step02_receptor_pairs_basin.m
%                                            rows = pairs, cols = node indices
%   INPUT/<...>.mat                         - urban dispersion data
%
% OUTPUTS:
%   M_triplets_<city>_<dir>.txt            - triplet metric matrix sorted by
%                                            R_abc (ascending)
%                                            cols: [a, b, c, cardinality, R_abc]
%   best_triplet_basin_<city>_<dir>_delta0<d>.txt  - sub-basins of the best
%                                            triplet (3 cols: ab, bc, ac),
%                                            NaN-padded to equal length
%   triplets_single_node_<city>_<dir>_delta0<d>.txt - individual node score
%                                            R_a summed over all triplets
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% USER PARAMETERS
% -------------------------------------------------------------------------

city_names   = {'Firenze', 'Parigi', 'New_York', 'Lione'};
phi_rotation = [34, 23, 124, 169];    % rotation angle of main street axis [deg]
wind_dirs    = [360, 45, 90, 135, 180, 225, 270, 315];  % wind directions [deg]
n_nodes_list = [706, 1198, 2460, 748];
delta        = 0.1;                    % relative measurement error [-]

addpath('FUNCTIONS')

% -------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------

for idx_city = 1 %choose the city in city_names

    n_nodes     = n_nodes_list(idx_city);

    for idx_dir = 7  %choose the wind direction in wind_dirs

        % Effective wind direction (accounting for street-axis rotation)
        wind_dir_effective = wind_dirs(idx_dir) + phi_rotation(idx_city);
        if wind_dir_effective > 360
            wind_dir_effective = wind_dir_effective - 360;
        end

        % -----------------------------------------------------------------
        % LOAD PAIR-LEVEL RESULTS (outputs of previous pipeline steps)
        % -----------------------------------------------------------------

        % Unsorted pair metric matrix: cols = [a, b, B_(ab,n), P_ab, R_ab]
        pair_metrics = table2array(readtable( ...
            ['M_' city_names{idx_city} '_' num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '_original.txt']));

        % Detection basin node lists: row k = basin nodes of pair k (zero-padded)
        pair_basins = table2array(readtable( ...
            ['pairs_' city_names{idx_city} '_' num2str(wind_dir_effective) '_basin_nodes.txt']));


        % -----------------------------------------------------------------
        % INITIALISE TRIPLET METRIC MATRIX
        % -----------------------------------------------------------------
        n               = n_nodes;
        n_triplets      = n * (n-1) * (n-2) / 6;
        M               = zeros(n_triplets, 5);
        counter_triplet = 0;

        % Precompute cumulative pair-row offsets for indexing into
        % pair_metrics and pair_basins.
        % y(d) = n - d is the number of pairs (d, b) with b > d.
        % cumsum_y(k) gives the total number of pairs with first index <= k,
        % i.e. the row offset needed to locate pairs starting at index k+1.
        y        = zeros(1, n-2);
        for d = 1 : n-2
            y(d) = n - d;
        end

        % =================================================================
        % TRIPLET LOOP — special case a = 1
        %
        % When a = 1, pair row indices simplify to:
        %   row(1, b) = b - 1
        %   row(1, c) = c - 1
        %   row(b, c) = cumsum(y(1:b-1)) + (c - b)
        % =================================================================
        for a = 1
            for b = a : n
                for c = b : n
                    if a ~= b && b ~= c && a ~= c

                        counter_triplet = counter_triplet + 1;

                        y_b = sum(y(1 : b-1));

                        % R values of the three constituent pairs
                        R_vec = [pair_metrics(b-a,         5), ...
                                 pair_metrics(c-a,         5), ...
                                 pair_metrics(y_b + (c-b), 5)];
                        [~, I_best_pair] = max(R_vec);

                        % Basin node lists for the three pairs
                        v_basin_ab = nonzeros(pair_basins(b-a,         :));
                        v_basin_ac = nonzeros(pair_basins(c-a,         :));
                        v_basin_bc = nonzeros(pair_basins(y_b + (c-b), :));

                        % Best pair keeps its full basin. Overlapping nodes
                        % are removed from the other two pairs' basins.
                        % A secondary comparison then resolves the residual
                        % overlap between those two pairs.
                        if I_best_pair == 1
                            % ab is dominant
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            v_basin_ac  = setdiff(v_basin_ac, inter_ab_ac);
                            inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                            v_basin_bc  = setdiff(v_basin_bc, inter_ab_bc);
                            [~, I_sec]  = max([pair_metrics(c-a, 5), pair_metrics(y_b+(c-b), 5)]);
                            inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                            if I_sec == 1
                                v_basin_bc = setdiff(v_basin_bc, inter_ac_bc);
                            elseif I_sec == 2
                                v_basin_ac = setdiff(v_basin_ac, inter_ac_bc);
                            end

                        elseif I_best_pair == 2
                            % ac is dominant
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            v_basin_ab  = setdiff(v_basin_ab, inter_ab_ac);
                            inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                            v_basin_bc  = setdiff(v_basin_bc, inter_ac_bc);
                            [~, I_sec]  = max([pair_metrics(b-a, 5), pair_metrics(y_b+(c-b), 5)]);
                            inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                            if I_sec == 1
                                v_basin_bc = setdiff(v_basin_bc, inter_ab_bc);
                            elseif I_sec == 2
                                v_basin_ab = setdiff(v_basin_ab, inter_ab_bc);
                            end

                        elseif I_best_pair == 3
                            % bc is dominant
                            inter_bc_ac = intersect(v_basin_bc, v_basin_ac);
                            v_basin_ac  = setdiff(v_basin_ac, inter_bc_ac);
                            inter_bc_ab = intersect(v_basin_ab, v_basin_bc);
                            v_basin_ab  = setdiff(v_basin_ab, inter_bc_ab);
                            [~, I_sec]  = max([pair_metrics(b-a, 5), pair_metrics(c-a, 5)]);
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            if I_sec == 1
                                v_basin_ac = setdiff(v_basin_ac, inter_ab_ac);
                            elseif I_sec == 2
                                v_basin_ab = setdiff(v_basin_ab, inter_ab_ac);
                            end
                        end

                        cardinality = length(v_basin_ab) + length(v_basin_ac) + length(v_basin_bc);
                        M(counter_triplet, :) = [a, b, c, cardinality, ...
                            pair_metrics(b-a, 5) + pair_metrics(c-a, 5) + pair_metrics(y_b+(c-b), 5)];

                    end
                end
            end
        end

        % =================================================================
        % TRIPLET LOOP — general case a = 2 .. n-1
        %
        % When a > 1, pair row indices require the cumulative offset y_a:
        %   row(a, b) = cumsum(y(1:a-1)) + (b - a)
        %   row(a, c) = cumsum(y(1:a-1)) + (c - a)
        %   row(b, c) = cumsum(y(1:b-1)) + (c - b)
        % =================================================================
        for a = 2 : (n-1)
            for b = a : n
                for c = b : n
                    if a ~= b && b ~= c && a ~= c

                        counter_triplet = counter_triplet + 1;

                        if mod(counter_triplet, 10000) == 0
                            counter_triplet/n_triplets*100
                        end

                        y_a = sum(y(1 : a-1));
                        y_b = sum(y(1 : b-1));

                        R_vec = [pair_metrics(y_a + (b-a), 5), ...
                                 pair_metrics(y_a + (c-a), 5), ...
                                 pair_metrics(y_b + (c-b), 5)];

                        [~, I_best_pair] = max(R_vec);

                        v_basin_ab = nonzeros(pair_basins(y_a + (b-a), :));
                        v_basin_ac = nonzeros(pair_basins(y_a + (c-a), :));
                        v_basin_bc = nonzeros(pair_basins(y_b + (c-b), :));

                        if I_best_pair == 1
                            % ab is dominant
                            % Note: secondary comparison uses column 3 (B index)
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            v_basin_ac  = setdiff(v_basin_ac, inter_ab_ac);
                            inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                            v_basin_bc  = setdiff(v_basin_bc, inter_ab_bc);
                            [~, I_sec]  = max([pair_metrics(y_a+(c-a), 3), pair_metrics(y_b+(c-b), 3)]);
                            inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                            if I_sec == 1
                                v_basin_bc = setdiff(v_basin_bc, inter_ac_bc);
                            elseif I_sec == 2
                                v_basin_ac = setdiff(v_basin_ac, inter_ac_bc);
                            end

                        elseif I_best_pair == 2
                            % ac is dominant
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            v_basin_ab  = setdiff(v_basin_ab, inter_ab_ac);
                            inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                            v_basin_bc  = setdiff(v_basin_bc, inter_ac_bc);
                            [~, I_sec]  = max([pair_metrics(y_a+(b-a), 5), pair_metrics(y_b+(c-b), 5)]);
                            inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                            if I_sec == 1
                                v_basin_bc = setdiff(v_basin_bc, inter_ab_bc);
                            elseif I_sec == 2
                                v_basin_ab = setdiff(v_basin_ab, inter_ab_bc);
                            end

                        elseif I_best_pair == 3
                            % bc is dominant
                            inter_bc_ac = intersect(v_basin_bc, v_basin_ac);
                            v_basin_ac  = setdiff(v_basin_ac, inter_bc_ac);
                            inter_bc_ab = intersect(v_basin_ab, v_basin_bc);
                            v_basin_ab  = setdiff(v_basin_ab, inter_bc_ab);
                            [~, I_sec]  = max([pair_metrics(y_a+(b-a), 5), pair_metrics(y_a+(c-a), 5)]);
                            inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                            if I_sec == 1
                                v_basin_ac = setdiff(v_basin_ac, inter_ab_ac);
                            elseif I_sec == 2
                                v_basin_ab = setdiff(v_basin_ab, inter_ab_ac);
                            end
                        end

                        cardinality = length(v_basin_ab) + length(v_basin_ac) + length(v_basin_bc);
                        M(counter_triplet, :) = [a, b, c, cardinality, ...
                            pair_metrics(y_a+(b-a), 5) + pair_metrics(y_a+(c-a), 5) + pair_metrics(y_b+(c-b), 5)];

                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % SORT TRIPLET METRIC MATRIX BY R_abc (column 5, ascending)
        % -----------------------------------------------------------------
        M_col1 = M(:, 1);
        M_col2 = M(:, 2);
        M_col3 = M(:, 3);
        M_col4 = M(:, 4);

        [R_sorted, I_sorted] = sort(M(:, 5));
        M(:, 1) = M_col1(I_sorted);
        M(:, 2) = M_col2(I_sorted);
        M(:, 3) = M_col3(I_sorted);
        M(:, 4) = M_col4(I_sorted);
        M(:, 5) = R_sorted;

        writematrix(M, ['M_triplets_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

        % =================================================================
        % CHARACTERISE THE BEST TRIPLET
        % =================================================================

        a_best = M(end, 1);
        b_best = M(end, 2);
        c_best = M(end, 3);

        if a_best == 1

            y_b = sum(y(1 : b_best-1));

            R_vec = [pair_metrics(b_best - a_best,         5), ...
                     pair_metrics(c_best - a_best,         5), ...
                     pair_metrics(y_b + (c_best - b_best), 5)];
            [~, I_best_pair] = max(R_vec);

            v_basin_ab = nonzeros(pair_basins(b_best - a_best,         :));
            v_basin_ac = nonzeros(pair_basins(c_best - a_best,         :));
            v_basin_bc = nonzeros(pair_basins(y_b + (c_best - b_best), :));

            if I_best_pair == 1
                inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                v_basin_ac  = setdiff(v_basin_ac, inter_ab_ac);
                inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                v_basin_bc  = setdiff(v_basin_bc, inter_ab_bc);
                [~, I_sec]  = max([pair_metrics(c_best-a_best, 5), pair_metrics(y_b+(c_best-b_best), 5)]);
                inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                if I_sec == 1
                    v_basin_bc = setdiff(v_basin_bc, inter_ac_bc);
                elseif I_sec == 2
                    v_basin_ac = setdiff(v_basin_ac, inter_ac_bc);
                end

            elseif I_best_pair == 2
                inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                v_basin_ab  = setdiff(v_basin_ab, inter_ab_ac);
                inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                v_basin_bc  = setdiff(v_basin_bc, inter_ac_bc);
                [~, I_sec]  = max([pair_metrics(b_best-a_best, 5), pair_metrics(y_b+(c_best-b_best), 5)]);
                inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                if I_sec == 1
                    v_basin_bc = setdiff(v_basin_bc, inter_ab_bc);
                elseif I_sec == 2
                    v_basin_ab = setdiff(v_basin_ab, inter_ab_bc);
                end

            elseif I_best_pair == 3
                inter_bc_ac = intersect(v_basin_bc, v_basin_ac);
                v_basin_ac  = setdiff(v_basin_ac, inter_bc_ac);
                inter_bc_ab = intersect(v_basin_ab, v_basin_bc);
                v_basin_ab  = setdiff(v_basin_ab, inter_bc_ab);
                [~, I_sec]  = max([pair_metrics(b_best-a_best, 5), pair_metrics(c_best-a_best, 5)]);
                inter_ab_ac = intersect(v_basin_ab, v_basin_bc);   % as in original
                if I_sec == 1
                    v_basin_ac = setdiff(v_basin_ac, inter_ab_ac);
                elseif I_sec == 2
                    v_basin_ab = setdiff(v_basin_ab, inter_ab_ac);
                end
            end

            P = [pair_metrics(b_best-a_best,         5) / length(v_basin_ab), ...
                 pair_metrics(c_best-a_best,         5) / length(v_basin_ac), ...
                 pair_metrics(y_b+(c_best-b_best),   5) / length(v_basin_bc)];

        else  % a_best > 1

            y_a = sum(y(1 : a_best-1));
            y_b = sum(y(1 : b_best-1));

            R_vec = [pair_metrics(y_a + (b_best-a_best), 5), ...
                     pair_metrics(y_a + (c_best-a_best), 5), ...
                     pair_metrics(y_b + (c_best-b_best), 5)];
            [~, I_best_pair] = max(R_vec);

            v_basin_ab = nonzeros(pair_basins(y_a + (b_best-a_best), :));
            v_basin_ac = nonzeros(pair_basins(y_a + (c_best-a_best), :));
            v_basin_bc = nonzeros(pair_basins(y_b + (c_best-b_best), :));

            if I_best_pair == 1
                % Intersections on raw (zero-padded) rows, as in original
                inter_ab_ac = intersect(pair_basins(y_a+(b_best-a_best), :), pair_basins(y_a+(c_best-a_best), :));
                v_basin_ac  = setdiff(v_basin_ac, inter_ab_ac);
                inter_ab_bc = intersect(pair_basins(y_a+(b_best-a_best), :), pair_basins(y_b+(c_best-b_best), :));
                v_basin_bc  = setdiff(v_basin_bc, inter_ab_bc);
                [~, I_sec]  = max([pair_metrics(y_a+(c_best-a_best), 3), pair_metrics(y_b+(c_best-b_best), 3)]);
                inter_ac_bc = intersect(pair_basins(y_a+(c_best-a_best), :), pair_basins(y_b+(c_best-b_best), :));
                if I_sec == 1
                    v_basin_bc = setdiff(v_basin_bc, inter_ac_bc);
                elseif I_sec == 2
                    v_basin_ac = setdiff(v_basin_ac, inter_ac_bc);
                end

            elseif I_best_pair == 2
                inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                v_basin_ab  = setdiff(v_basin_ab, inter_ab_ac);
                inter_ac_bc = intersect(v_basin_ac, v_basin_bc);
                v_basin_bc  = setdiff(v_basin_bc, inter_ac_bc);
                [~, I_sec]  = max([pair_metrics(y_a+(b_best-a_best), 5), pair_metrics(y_b+(c_best-b_best), 5)]);
                inter_ab_bc = intersect(v_basin_ab, v_basin_bc);
                if I_sec == 1
                    v_basin_bc = setdiff(v_basin_bc, inter_ab_bc);
                elseif I_sec == 2
                    v_basin_ab = setdiff(v_basin_ab, inter_ab_bc);
                end

            elseif I_best_pair == 3
                inter_bc_ac = intersect(v_basin_bc, v_basin_ac);
                v_basin_ac  = setdiff(v_basin_ac, inter_bc_ac);
                inter_bc_ab = intersect(v_basin_ab, v_basin_bc);
                v_basin_ab  = setdiff(v_basin_ab, inter_bc_ab);
                [~, I_sec]  = max([pair_metrics(y_a+(b_best-a_best), 5), pair_metrics(y_a+(c_best-a_best), 5)]);
                inter_ab_ac = intersect(v_basin_ab, v_basin_ac);
                if I_sec == 1
                    v_basin_ac = setdiff(v_basin_ac, inter_ab_ac);
                elseif I_sec == 2
                    v_basin_ab = setdiff(v_basin_ab, inter_ab_ac);
                end
            end

            P = [pair_metrics(y_a+(b_best-a_best), 5) / length(v_basin_ab), ...
                 pair_metrics(y_a+(c_best-a_best), 5) / length(v_basin_ac), ...
                 pair_metrics(y_b+(c_best-b_best), 5) / length(v_basin_bc)];

        end

        % Write detection sub-basins as three columns, NaN-padded to equal length
        v_cell  = {v_basin_ab, v_basin_bc, v_basin_ac};
        max_len = max(cellfun(@numel, v_cell));
        v_basins = cell2mat(cellfun( ...
            @(v) [v(:); nan(max_len - numel(v), 1)], v_cell, 'UniformOutput', false));

        writematrix(v_basins, ['best_triplet_basin_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

        % =================================================================
        % INDIVIDUAL NODE PERFORMANCE SCORE
        % R_a = sum of R_abc over all triplets that include node a
        % =================================================================

        M_nodes = zeros(n, 2);

        for g = 1 : n
            S_g      = 0;
            idx_as_a = find(M(:, 1) == g);
            idx_as_b = find(M(:, 2) == g);
            idx_as_c = find(M(:, 3) == g);

            for f = 1 : length(idx_as_a)
                S_g = S_g + M(idx_as_a(f), 5);
            end
            for f = 1 : length(idx_as_b)
                S_g = S_g + M(idx_as_b(f), 5);
            end
            for f = 1 : length(idx_as_c)
                S_g = S_g + M(idx_as_c(f), 5);
            end

            M_nodes(g, :) = [g, S_g];
        end

        writematrix(M_nodes, ['triplets_single_node_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

    end % wind direction loop
end     % city loop