% =========================================================================
% RECEPTOR PAIR IDENTIFICATION FOR URBAN POLLUTANT SOURCE DETECTION
% =========================================================================
%
% This script identifies optimal pairs of receptor nodes in urban road
% networks for pollutant source identification, based on the matrix method
% described in:
%
%   Carlini, B., Salizzoni, P., Ridolfi, L., Fellini, S. (2026).
%   A model-independent matrix method for optimizing urban air quality
%   monitoring networks. Environmental Modelling & Software.
%
% For each city and wind direction, it computes:
%
%   - The detection basin index  B_(ab,n) : fraction of nodes detectable
%                                           by the receptor pair (a,b)
%   - The precision index        P_ab     : discriminability between source
%                                           nodes within the detection basin
%   - The combined performance   R_ab     : overall quality of pair (a,b)
%
% The best receptor pair is identified, and its detection basin is
% partitioned into certain nodes (nodes uniquely localised as sources) and uncertain nodes
% (nodes partially discriminated as sources).
%
% Additionally, a scalar performance score R_a is computed for each
% individual node, summarising its contribution across all pairs.
%
% -------------------------------------------------------------------------
% NOTE ON MODEL INDEPENDENCE
%
% The sensor placement algorithm operates exclusively on the transport
% matrix A (see Section 2 of the paper), which encodes the source-receptor
% relationships across the urban domain. In this script, A is constructed
% using the simplified network-based dispersion model of Fellini et al.
% (2019, 2020) as an illustrative example (see Appendix A of the paper).
%
% The construction of A (see Section "BUILD TRANSPORT MATRIX A" below) is
% entirely decoupled from the sensor placement algorithm. A can be replaced
% by any source-receptor matrix derived from alternative dispersion models
% (e.g., CFD simulations, Gaussian models) or from observational data,
% without any modification to the algorithm.
% -------------------------------------------------------------------------
%
% INPUTS  (loaded from INPUT\ folder as .mat files):
%   G_str          - directed graph of the urban road network (MATLAB graph)
%   L_str_ord_a    - link lengths [m]
%   Ud_str_ord_a   - mean wind speed projected on each link [m/s]
%   H_str_ord_a    - mean building height along each link [m]
%   U_str_ord_a    - mean street-canyon wind speed [m/s]
%
% OUTPUTS (written as .txt files in the working directory):
%   M_<city>_<dir>_delta0<d>_original.txt       - full metric matrix (unsorted)
%   M_bestpair_<city>_<dir>_delta0<d>.txt        - metric matrix sorted by R_ab
%   v_certain_nodes_<city>_<dir>_delta0<d>.txt   - certain source nodes (Pi = 1)
%   v_uncertain_<city>_<dir>_delta0<d>.txt       - uncertain source nodes (Pi > 0)
%   v_zero_nodes_<city>_<dir>_delta0<d>.txt      - basin nodes with Pi = 0
%   single_node_<city>_<dir>_delta0<d>.txt       - individual node performance
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% USER PARAMETERS
% -------------------------------------------------------------------------

city_names   = {'Firenze', 'Parigi', 'New_York', 'Lione'};  % city labels
phi_rotation = [34, 23, 124, 169];    % rotation angle of main street axis [deg]
wind_dirs    = [360, 45, 90, 135, 180, 225, 270, 315];  % wind directions [deg]
n_nodes_list = [706, 1198, 2460, 748]; % number of nodes per city
delta        = 0.1;                    % relative measurement error [-]

addpath('FUNCTIONS')

% -------------------------------------------------------------------------
% MAIN LOOP OVER CITIES AND WIND DIRECTIONS
% -------------------------------------------------------------------------

for idx_city = 1 %choose the city in city_names

    n_nodes     = n_nodes_list(idx_city);
    U_str_ord_a = zeros(n_nodes, 1);

    for idx_dir = 7 %choose the wind direction in wind_dirs

        % Compute effective wind direction (accounting for street-axis rotation)
        wind_dir_effective = wind_dirs(idx_dir) + phi_rotation(idx_city);
        if wind_dir_effective > 360
            wind_dir_effective = wind_dir_effective - 360;
        end

        % -----------------------------------------------------------------
        % LOAD INPUT DATA FOR CITY GRAPH 
        % -----------------------------------------------------------------
        input_file = fullfile('INPUT', ...
            ['5_', num2str(wind_dir_effective), '_C0_Cth_10_', ...
             city_names{idx_city}, '_25_03.mat']);
        load(input_file)

        % -----------------------------------------------------------------
        % BUILD TRANSPORT MATRIX A
        %
        % This section constructs the source-receptor transport matrix A
        % using the network-based dispersion model of Fellini et al. (2019,
        % 2020), described in Appendix A of the paper. Each element A(i,j)
        % represents the normalised concentration at receptor i due to a
        % unit emission at source j, computed as the product of exponential
        % decay weights along the shortest path between j and i on the
        % urban street network.
        %
        % To use an alternative dispersion model, replace this section with
        % any routine that produces an (n x n) matrix A where A(i,j) > 0
        % if source j influences receptor i, and A(i,j) = 0 otherwise.
        % The sensor placement algorithm below requires no further changes.
        % -----------------------------------------------------------------

        n_links = length(U_str_ord_a);

        % Replace zero wind speeds with a small value to avoid division by zero
        U_str_ord_a(U_str_ord_a == 0) = 1e-4;

        % Dimensionless dispersion weight for each link (Eq. A.3 in the paper)
        w_link = L_str_ord_a .* Ud_str_ord_a ./ (H_str_ord_a .* U_str_ord_a);

        % Compute exponential decay matrix via shortest-path algorithm
        % (distances_threshold.m, which calls dijkstra_threshold.m):
        % D_decay(i,j) = product of decay weights along shortest path from i to j
        % Transposing gives A(i,j) = influence of source j on receptor i
        D_decay = distances_threshold(G_str, w_link);
        A = D_decay.';
        A(A < 1e-6) = 0;   % threshold small values to enforce sparsity

        % -----------------------------------------------------------------
        % INITIALISE METRIC MATRIX M
        %
        % Columns of M:
        %   1 - receptor node index a
        %   2 - receptor node index b
        %   3 - normalised detection basin  B_(ab,n) = cardinality / n
        %   4 - average precision index     P_ab
        %   5 - combined performance index  R_ab
        % -----------------------------------------------------------------

        n           = n_nodes;
        n_pairs     = n * (n - 1) / 2;
        M           = zeros(n_pairs, 5);
        counter_pair = 0;

        % -----------------------------------------------------------------
        % COMPUTE METRICS FOR ALL RECEPTOR PAIRS (a, b)
        % -----------------------------------------------------------------

        for a = 1 : (n - 1)
            for b = (a + 1) : n

                counter_pair = counter_pair + 1;

                % Ratio vectors (indexed by source node j)
                ratios = zeros(1, n);
                left   = zeros(1, n);
                right  = zeros(1, n);
                cardinality = 0;

                % --- Detection basin: nodes j reachable by both a and b ---
                for j = 1 : n
                    if A(a, j) ~= 0 && A(b, j) ~= 0
                        cardinality  = cardinality + 1;
                        ratios(j)    = A(a, j) / A(b, j);
                        left(j)      = ((1 - delta) * A(a, j)) / ((1 + delta) * A(b, j));
                        right(j)     = ((1 + delta) * A(a, j)) / ((1 - delta) * A(b, j));
                    end
                end

                % Sort ratios to obtain ordered index i (from j)
                % (see Section 2.3 and Eq. 7 of the paper)
                [ratios_sorted, ord_i] = sort(ratios);
                left_sorted  = left(ord_i);
                right_sorted = right(ord_i);

                ratios_nonzero = nonzeros(ratios_sorted);
                metric_free_area = 0;

                % --- Handle degenerate cases ---
                if isempty(ratios_nonzero)
                    % No common detectable nodes
                    M(counter_pair, :) = [a, b, 0, 0, 0];

                elseif length(ratios_nonzero) > 1 && all(ratios_nonzero(1) == ratios_nonzero(2:end))
                    % All ratios identical: detection possible but no discrimination
                    M(counter_pair, :) = [a, b, cardinality / n, 0, 0];

                else
                    % --- Precision index computation ---
                    % Computes the reliable interval RI_i for each source node i
                    % (Eq. 9 of the paper) and accumulates the precision metric P_ab
                    % (Eq. 12-13 of the paper)
                    left_nonzero  = nonzeros(left_sorted);
                    right_nonzero = nonzeros(right_sorted);

                    % Build signed flag array for interval tracking
                    counter_flag = 1 : n;
                    N_flag = left_nonzero;
                    N_flag(end+1 : end*2) = nonzeros(right_sorted);
                    idx_zero_limit = find(left_sorted == 0);
                    flag_a = (-1) * counter_flag(idx_zero_limit(end)+1 : n);
                    flag_a(end+1 : end*2) = -flag_a;
                    [~, I_sort_flag] = sort(N_flag);
                    flag_a = flag_a(I_sort_flag);

                    for i = 1 : length(ratios_nonzero)
                        free_area = 0;

                        % Check if interval i is not dominated by any other interval
                        u = find(left_nonzero(:) <= left_nonzero(i));
                        u(u == i) = [];

                        if sum(right_nonzero(u) >= right_nonzero(i)) == 0

                            % Extract flag segment for node i
                            flag_e = flag_a( ...
                                find(flag_a == -(i + idx_zero_limit(end))) : ...
                                find(flag_a ==  (i + idx_zero_limit(end))));
                            B_a = flag_e;

                            % Prune flag array: forward and backward passes
                            for e = 1 : (length(flag_e) - 2)

                                % Forward pass
                                if flag_e(e+1) > 0
                                    x = find(flag_e == -flag_e(e+1));
                                    if ~isempty(x)
                                        if e >= x + 1
                                            B_a(x+1 : e) = 0;
                                        end
                                    else
                                        B_a(1 : e) = 0;
                                    end
                                end

                                % Backward pass
                                if flag_e(end-e) < 0
                                    y = find(flag_e == -flag_e(end-e));
                                    if ~isempty(y)
                                        if (y-1) >= (length(flag_e) - e + 1)
                                            B_a(end-e+1 : y-1) = 0;
                                        end
                                    else
                                        B_a(end-e+1 : end) = 0;
                                    end
                                end

                            end

                            % Accumulate free area from non-overlapping interval segments
                            B_a_nonzero = nonzeros(B_a);
                            for k = 1 : 2 : length(B_a_nonzero)
                                bk  = B_a_nonzero(k);
                                bk1 = B_a_nonzero(k+1);
                                if     bk  < 0 && bk1 < 0
                                    free_area = free_area + (left_sorted(-bk1)  - left_sorted(-bk));
                                elseif bk  > 0 && bk1 > 0
                                    free_area = free_area + (right_sorted(bk1)  - right_sorted(bk));
                                elseif bk  > 0 && bk1 < 0
                                    free_area = free_area + (left_sorted(-bk1)  - right_sorted(bk));
                                elseif bk  < 0 && bk1 > 0
                                    free_area = free_area + (right_sorted(bk1)  - left_sorted(-bk));
                                end
                            end

                            metric_free_area = metric_free_area + ...
                                free_area / (right_nonzero(i) - left_nonzero(i));

                        end
                    end

                    % Store metrics: B_(ab,n), P_ab, R_ab (Eqs. 11, 13, 14)
                    M(counter_pair, :) = [a, b, ...
                        cardinality / n, ...
                        metric_free_area / cardinality, ...
                        metric_free_area / n];

                end
            end
        end

        % Save original (unsorted) metric matrix
        writematrix(M, ['M_' city_names{idx_city} '_' num2str(wind_dir_effective) ...
            '_delta0' num2str(delta * 10) '_original.txt'])

        % -----------------------------------------------------------------
        % SORT METRIC MATRIX BY COMBINED PERFORMANCE INDEX R_ab (column 5)
        % -----------------------------------------------------------------

        % Keep original column copies before sorting
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

        writematrix(M, ['M_bestpair_' city_names{idx_city} '_' num2str(wind_dir_effective) ...
            '_delta0' num2str(delta * 10) '.txt'])

        % =================================================================
        % CHARACTERISE THE BEST RECEPTOR PAIR
        %
        % For the pair (a*, b*) with highest R_ab, classify each node in
        % the detection basin as:
        %   - certain   : precision index = 1  (unique source localisation)
        %   - uncertain : precision index in (0, 1)  (partial discrimination)
        %   - Pi = 0    : not included in either output vector
        % =================================================================

        a_best = M(end, 1);
        b_best = M(end, 2);

        ratios = zeros(1, n);
        left   = zeros(1, n);
        right  = zeros(1, n);
        cardinality = 0;
        basin_nodes = [];   % nodes in the detection basin of the best pair

        for j = 1 : n
            if A(a_best, j) ~= 0 && A(b_best, j) ~= 0
                cardinality      = cardinality + 1;
                basin_nodes(end+1) = j;                             %#ok<AGROW>
                ratios(j)          = A(a_best, j) / A(b_best, j);
                left(j)            = ((1 - delta) * A(a_best, j)) / ((1 + delta) * A(b_best, j));
                right(j)           = ((1 + delta) * A(a_best, j)) / ((1 - delta) * A(b_best, j));
            end
        end

        [ratios_sorted, ord_i] = sort(ratios);
        left_sorted  = left(ord_i);
        right_sorted = right(ord_i);

        ratios_nonzero = nonzeros(ratios_sorted);
        metric_free_area = 0;

        certain_nodes   = zeros(1, length(ratios_nonzero));
        uncertain_nodes = zeros(1, length(ratios_nonzero));
        n_certain       = 0;
        n_uncertain     = 0;

        if ~isempty(ratios_nonzero) && ~(length(ratios_nonzero) > 1 && all(ratios_nonzero(1) == ratios_nonzero(2:end)))

            left_nonzero  = nonzeros(left_sorted);
            right_nonzero = nonzeros(right_sorted);

            counter_flag     = 1 : n;
            N_flag           = left_nonzero;
            N_flag(end+1 : end*2) = nonzeros(right_sorted);
            idx_zero_limit   = find(left_sorted == 0);
            flag_a           = (-1) * counter_flag(idx_zero_limit(end)+1 : n);
            flag_a(end+1 : end*2) = -flag_a;
            [~, I_sort_flag] = sort(N_flag);
            flag_a           = flag_a(I_sort_flag);

            for i = 1 : length(ratios_nonzero)
                free_area = 0;

                u = find(left_nonzero(:) <= left_nonzero(i));
                u(u == i) = [];

                if sum(right_nonzero(u) >= right_nonzero(i)) == 0

                    flag_e = flag_a( ...
                        find(flag_a == -(i + idx_zero_limit(end))) : ...
                        find(flag_a ==  (i + idx_zero_limit(end))));
                    B_a = flag_e;

                    for e = 1 : (length(flag_e) - 2)

                        % Forward pass
                        if flag_e(e+1) > 0
                            x = find(flag_e == -flag_e(e+1));
                            if ~isempty(x)
                                if e >= x + 1
                                    B_a(x+1 : e) = 0;
                                end
                            else
                                B_a(1 : e) = 0;
                            end
                        end

                        % Backward pass
                        if flag_e(end-e) < 0
                            y = find(flag_e == -flag_e(end-e));
                            if ~isempty(y)
                                if (y-1) >= (length(flag_e) - e + 1)
                                    B_a(end-e+1 : y-1) = 0;
                                end
                            else
                                B_a(end-e+1 : end) = 0;
                            end
                        end

                    end

                    B_a_nonzero = nonzeros(B_a);
                    for k = 1 : 2 : length(B_a_nonzero)
                        bk  = B_a_nonzero(k);
                        bk1 = B_a_nonzero(k+1);
                        if     bk  < 0 && bk1 < 0
                            free_area = free_area + (left_sorted(-bk1)  - left_sorted(-bk));
                        elseif bk  > 0 && bk1 > 0
                            free_area = free_area + (right_sorted(bk1)  - right_sorted(bk));
                        elseif bk  > 0 && bk1 < 0
                            free_area = free_area + (left_sorted(-bk1)  - right_sorted(bk));
                        elseif bk  < 0 && bk1 > 0
                            free_area = free_area + (right_sorted(bk1)  - left_sorted(-bk));
                        end
                    end

                    precision_i = free_area / (right_nonzero(i) - left_nonzero(i));

                    if precision_i == 1
                        % Source node i is uniquely identifiable
                        n_certain = n_certain + 1;
                        certain_nodes(n_certain) = ord_i(i + idx_zero_limit(end));
                    elseif precision_i > 0
                        % Source node i is partially discriminated
                        n_uncertain = n_uncertain + 1;
                        uncertain_nodes(n_uncertain) = ord_i(i + idx_zero_limit(end));
                    end

                end
            end

        end

        writematrix(certain_nodes,  ['v_certain_nodes_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])
        writematrix(uncertain_nodes, ['v_uncertain_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

        % Nodes in the detection basin with Pi = 0 (detectable but not
        % discriminable): these are basin nodes that appear in neither
        % certain_nodes nor uncertain_nodes. Used by ELABORATION/step04_plot_basin.m
        % to colour these nodes in blue (as in Figure 4 of the paper).
        certain_nodes_valid  = certain_nodes(certain_nodes > 0);
        uncertain_nodes_valid = uncertain_nodes(uncertain_nodes > 0);
        zero_nodes = basin_nodes(~ismember(basin_nodes, ...
            [certain_nodes_valid, uncertain_nodes_valid]));
        writematrix(zero_nodes, ['v_zero_nodes_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

        % =================================================================
        % IDENTIFY ALL BEST PAIRS (tied R_ab) AND INDIVIDUAL NODE SCORES
        % =================================================================

        % --- All pairs sharing the maximum R_ab value ---
        best_R        = M(end, 5);
        idx_best      = find(M(:, 5) == best_R);
        M_best_pairs  = zeros(length(idx_best), 2);

        for z = 1 : length(idx_best)
            M_best_pairs(z, :) = [M(end - z + 1, 1), M(end - z + 1, 2)];
        end

        % --- Individual node performance score R_a (Eq. 15 of the paper) ---
        % R_a is the sum of R_ab over all pairs (a, b) that include node a.
        % It summarises how useful node a is as a receptor across all pairings.

        M_nodes = zeros(n, 2);

        for g = 1 : n
            S_g = 0;

            % Pairs where g appears as the first receptor
            idx_as_a = find(M(:, 1) == g);
            S_g = S_g + sum(M(idx_as_a, 5));

            % Pairs where g appears as the second receptor
            idx_as_b = find(M(:, 2) == g);
            S_g = S_g + sum(M(idx_as_b, 5));

            M_nodes(g, :) = [g, S_g];
        end

        writematrix(M_nodes, ['single_node_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_delta0' num2str(delta * 10) '.txt'])

    end % wind direction loop
end     % city loop