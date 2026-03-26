% =========================================================================
% DETECTION BASIN NODE LISTS FOR RECEPTOR PAIRS
% =========================================================================
%
% This script is the second step in the receptor identification pipeline,
% and must be run after step01_receptor_pairs.m.
%
% For each receptor pair (a, b), it identifies the set of source nodes j
% that belong to the detection basin of the pair, i.e. nodes reachable by
% non-zero dispersion from both receptor a and receptor b.
%
% The results are stored in a matrix where each row corresponds to one
% receptor pair (ordered consistently with the metric matrix M produced by
% step01_receptor_pairs.m) and each column stores a basin node index
% (zero-padded to the maximum basin size across all pairs).
%
% This output is required as input by step03_receptor_triplets.m.
%
% -------------------------------------------------------------------------
% INPUTS  (loaded from INPUT\ folder as .mat files):
%   G_str          - directed graph of the urban road network
%   L_str_ord_a    - link lengths [m]
%   Ud_str_ord_a   - mean wind speed projected on each link [m/s]
%   H_str_ord_a    - mean building height along each link [m]
%   U_str_ord_a    - mean street-canyon wind speed [m/s]
%
% OUTPUTS:
%   pairs_<city>_<dir>_basin_nodes.txt  - matrix of basin node indices
%                                         rows  = receptor pairs (a, b)
%                                         cols  = node indices (zero-padded)
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% USER PARAMETERS
% -------------------------------------------------------------------------

city_names   = {'Firenze', 'Parigi', 'New_York', 'Lione'};
phi_rotation = [34, 23, 124, 169];   % rotation angle of main street axis [deg]
wind_dirs    = [360, 45, 90, 135, 180, 225, 270, 315];  % wind directions [deg]
n_nodes_list = [706, 1198, 2460, 748];

addpath('FUNCTIONS')

% -------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------

for idx_city = 1 %choose the city in city_names

    n_nodes     = n_nodes_list(idx_city);
    U_str_ord_a = zeros(n_nodes, 1);

    for idx_dir = 7 %choose the wind direction in wind_dirs

        % Effective wind direction
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
        % BUILD BASIN NODE MATRIX
        %
        % M_basins(pair_idx, :) contains the node indices j that belong to
        % the detection basin of pair (a, b), zero-padded to n columns.
        % Row ordering matches the metric matrix from
        % step01_receptor_pairs.m: pairs are enumerated as
        % (1,2),(1,3),...,(1,n),(2,3),...,(n-1,n).
        % -----------------------------------------------------------------
        n        = n_nodes;
        n_pairs  = n * (n-1) / 2;
        M_basins = zeros(n_pairs, n);
        counter_pair = 0;

        for a = 1 : (n-1)
            for b = (a+1) : n

                counter_pair = counter_pair + 1;
                counter_j    = 0;

                % A node j is in the detection basin of (a,b) if it has
                % non-zero dispersion connectivity to both a and b
                for j = 1 : n
                    if A(a, j) ~= 0 && A(b, j) ~= 0
                        counter_j = counter_j + 1;
                        M_basins(counter_pair, counter_j) = j;
                    end
                end

            end
        end

        % -----------------------------------------------------------------
        % WRITE OUTPUT
        % -----------------------------------------------------------------
        output_file = ['pairs_' city_names{idx_city} '_' ...
            num2str(wind_dir_effective) '_basin_nodes.txt'];
        writematrix(M_basins, output_file)

    end % wind direction loop
end     % city loop