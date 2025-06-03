display("In this file")
% Version 1: Initial version

% Version 2: - Add domains
%            - Add gaussian initial distribution of M-level population
%
% Version 3: Add separate population of m=+1/1
%
% Version 4: - Allow separate local dephasing rates for m=+1/-1 and m=0
%            - Allow separate cooperativities for m=+1/-1 and m=0
%            - refactor to handle case where number of collective atoms is zero
%
% Version 5: - If there's only one atom in the collective space for a domain
%              treat it as incoherent, not collective.
%            - Finally fix issues related to rounding domains to size of zero
% 
% Version 7: - Include ISC dark rates in the collective decay space by
%              adding a -(J+M)*gamma_isc to every population |J,M>
%
% Version 8: Use a more correct (i.e. sharper) detector response. The previous
%            version seemed to have characterised the IDQ using a pulse that
%            was too broad. FFS.


% matlab script to compare experimental NV fluorescence rates with
% simulation

% Assume that decoherence from single atom phase flips doesn't
% remove population from collective subspace and into non-collective
% subspace. Rather, it leaves it in the |J,M> state with probability
% |M/J|^2, and with probability 1-|M/J|^2 takes it into a seperate
% collective subspace with max spin J-1. I.e. it goes from |J,M>
% into |J-1, M-1>.

% Assume all population starts in the highest collective J subspace.
% Since we only have losses from this one, we can solve it for the
% entire time of interest. Then we can solve for the |J-1> subspace
% with losses into the |J-2> subspace, plus gains from losses coming
% from the |J> subspace for which we've already solved for the populations. 

% Also, for clarity, drop the possibility that we can be in the 
% +1 magnetic sublevel, and consider only the |e, m=0> -> |g, m=0> 
% manifold.

% Finally, drop the possibility of having multiple domains, as this
% will slow things down even more.

clear all
%close all


% Load the experimental data 

% Data is formatted as follows:
% There is data for two diamonds: NV1, 13, 14, 17, 36
% Each diamond has measurements at 8 different powers: 
% 1, 2, 5, 10, 20, 30, 40, 50uW

% Each measurement consists of the times (first column) followed by the counts (2nd column)
% So 80 columns ( 5 diamonds * 8 powers x 2 columns) in all

% Skip this since we're only doing theory curves
load 'cleaned_data_off_resonant.mat' 

num_datasets = size(fluorescence_counts_experiment, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation parameters
% Only these need to be changed



% For the paper we need to do NVs 1 (fitted), 17 (fitted), 29 (fitted), 56

% NV	N_max	gamma/2pi	gamma0_d/2pi	gamma1_d	Sigma0
% NV 01	  2	2.5 MHz		27 MHz		270 MHz		0.56
% NV 17   7	4.8 MHz		20 MHz		260 MHz		0.51
% NV 29   10	3.3 MHz		39 MHz		420 MHz		0.50
% NV 56	  50	7.9 MHz		20 MHz		450 MHz		0.50



diamond_number = 17;

% At most one of the three options below can be true
% If all are false we're on resonant and doing a partial Rabi flop for the excitation,
% with the degree of inversion set by proportion_excited_atoms

offresonant = 1;            % If true, the system is fully inverted
equal_populations = 0;      % If true, all the M-levels are equally populated
gaussian_populations = 0;   % If true, initial populations of the M-levels is gaussian centred about M=0;

gaussian_population_width = 0.4; % 1 std dev of total M-range, only used if gaussian_populations = 1
proportion_excited_atoms = 0.7;   % Only used if none of the above are true

N0 = 7; % Total number of atoms in m=0 magnetic sublevel
N1 = 0; % Total number of atoms in m=+1/-1 magnetic sublevel
 
num_domains = 1;
cooperativity_0 = 1;
cooperativity_1 = 1;

LDOSmultiplier = 3;
lifetime_zpl = LDOSmultiplier * 150e-9;     % assuming 5% branching ratio
lifetime_sideband = LDOSmultiplier * 12e-9;

lifetime_dark_0 = 88e-9; % Mott-seitz?
lifetime_dark_1 = 17e-9; % ISC rate, should be ~17ns at room temp according to literature

dephasing_lifetime_0 = 8e-9;
dephasing_lifetime_1 = 0.6e-9;

% We will allow a non-radiative decay acting on the upper m=0 state
% as per Monticone et al that gives fluorescence of the form
%
% I(t) = N gamma exp(-gamma t) exp(-k (gamma t)^c)
%k = 0;
%c = 0.5;   % 1/2 for dipole-dipole, 3/8 for quadrupole



non_rad_transition_is_collective = 0;

num_simulation_time_points = 7907;  %7907;
collective_effect_cutoff_time = 45e-9; % Collective part of simulation is really slow, so
                                       % stop doing after a certain time and assume all
                                       % remaining decays are non-collective

origin_shift = 25.36e-9;


% End simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Hack code for plotting several experimental traces
% slow: 1,3
% fast 2,8,17
% ultra fast 14, 36, 37, 56, 81, 97
%nv_index = 1;
%t_experiment_1 = 1e-9 * times_experiment(nv_index, :);
%fluorescence_counts_experiment_1 = fluorescence_counts_experiment(nv_index, :);
%peak_1 = max(fluorescence_counts_experiment_1);
%
%nv_index = 8;
%t_experiment_2 = 1e-9 * times_experiment(nv_index, :);
%fluorescence_counts_experiment_2 = fluorescence_counts_experiment(nv_index, :);
%peak_2 = max(fluorescence_counts_experiment_2);
%
%nv_index = 14;
%t_experiment_3 = 1e-9 * times_experiment(nv_index, :);
%fluorescence_counts_experiment_3 = fluorescence_counts_experiment(nv_index, :);
%peak_3 = max(fluorescence_counts_experiment_3);
%
%close all
%figure(1)
%hold on
%plot(t_experiment_1, fluorescence_counts_experiment_1 / peak_1, 'r', 'Linewidth', 2)
%plot(t_experiment_2, fluorescence_counts_experiment_2 / peak_2, 'g', 'Linewidth', 2)
%plot(t_experiment_3, fluorescence_counts_experiment_3 / peak_3, 'b', 'Linewidth', 2) 
%hold off


% Get the time data for this diamond and zero it sensibly
t_experiment = 1e-9 * times_experiment(diamond_number, :);
t_experiment = t_experiment - t_experiment(1);
number_of_time_bins = length(t_experiment);

start_time = t_experiment(1);
end_time = t_experiment(number_of_time_bins);
simulation_times = start_time : (end_time-start_time) / (num_simulation_time_points-1) : end_time;
dt = simulation_times(2) - simulation_times(1);


% We do the entire simulation twice, once for m=0 and once for m=+1/-1
% The variables total_* refer to the sum of sublevels 0 and +1/-1
% The variables domain_* refer to the per-domain quantities

total_visible_scattering_rate = zeros(1, num_simulation_time_points);
total_visible_collective_scattering_rate = zeros(1, num_simulation_time_points);
total_visible_noncollective_scattering_rate = zeros(1, num_simulation_time_points);

for sublevel = 0:1

  fprintf('\n\nDoing sublevel %i\n', sublevel);

  visible_scattering_rate = zeros(1, num_simulation_time_points);
  visible_collective_scattering_rate = zeros(1, num_simulation_time_points);
  visible_noncollective_scattering_rate = zeros(1, num_simulation_time_points);

  if sublevel == 0
    N = N0;
    cooperativity = cooperativity_0;
  else
    N = N1;
    cooperativity = cooperativity_1;
  end

  % if one of the sublevels has zero centres, just skip it
  if N==0
    continue;
  end


  %initial_collective_population = max(1, ceil(cooperativity * N));
  %initial_noncollective_population = N - initial_collective_population;

  % Apportion the centres over the domains. If, say, 5% of the cooperatively
  % acting centres are in domain 1, put 5% of the non-cooperatively acting
  % ones there too.

  clear centres_per_domains

  % Flat distribution of centres, i.e. same number of centres per domain
  for n=1:num_domains
    centres_per_domain(n) = ceil(N/num_domains);
  end

  %centres_per_domain
  % Gaussian distribution of number of centres per domain.
  %sigma = num_domains/4;
  %middle = (num_domains+1)/2;
  %for n=1:num_domains
  %  centres_per_domain(n) = round(N/(sigma*sqrt(2*pi)) * exp(-(n-middle)^2 / (2*sigma^2)));
  %end
  
%  if sublevel == 0
%    centres_per_domain(1) = 50;
%    centres_per_domain(2) = 1;
%    centres_per_domain(3) = 2;
%    centres_per_domain(4) = 2;
%    centres_per_domain(5) = 1;
%    centres_per_domain(6) = 1;
%    centres_per_domain(7) = 1;
%    centres_per_domain(8) = 10;
%    centres_per_domain(9) = 10;
%    centres_per_domain(10) = 10;
%    centres_per_domain(11) = 10;
%    centres_per_domain(12) = 10;
%    centres_per_domain(13) = 10;
%    centres_per_domain(14) = 10;
%    centres_per_domain(15) = 10;
%    centres_per_domain(16) = 15;
%    centres_per_domain(17) = 15;
%    centres_per_domain(18) = 45;
%  else
%    centres_per_domain(1) = 50;
%    centres_per_domain(2) = 1;
%    centres_per_domain(3) = 1;
%    centres_per_domain(4) = 1;
%    centres_per_domain(5) = 1;
%    centres_per_domain(6) = 1;
%    centres_per_domain(7) = 1;
%    centres_per_domain(8) = 10;
%    centres_per_domain(9) = 10;
%    centres_per_domain(10) = 10;
%    centres_per_domain(11) = 10;
%    centres_per_domain(12) = 10;
%    centres_per_domain(13) = 10;
%    centres_per_domain(14) = 10;
%    centres_per_domain(15) = 10;
%    centres_per_domain(16) = 15;
%    centres_per_domain(17) = 15;
%    centres_per_domain(18) = 45;
%  end

  %centres_per_domain(6) = N * 0.1;
  %centres_per_domain(7) = 0;
  %centres_per_domain(6) = 0;
  %centres_per_domain(7) = 0;

  % In each domain, assign collective and non-collective populations
  % Note these populations refer to the excited state only, not the ground state
  % We don't weight the number in the collective state by the proportion that
  % are excited here, because we take care of that later
  clear initial_collective_populations initial_noncollective_populations

  for n=1:num_domains
    initial_collective_populations(n) = floor(centres_per_domain(n)*cooperativity);
    initial_noncollective_populations(n) = floor(proportion_excited_atoms*(centres_per_domain(n) - initial_collective_populations(n)));
  end


  gamma_zpl = 1 / lifetime_zpl;   % decay rate from upper state into the ZPL
  gamma_sideband = 1 / lifetime_sideband;   %decay rate from upper state into phonon sidebands

  if sublevel == 0
    gamma_dark = 1/lifetime_dark_0;
    gamma = gamma_zpl + gamma_sideband + gamma_dark;   % total decay rate from upper state m=0
    gamma_dephasing = 1 / dephasing_lifetime_0;   % decay rate due to local dephasing for m=0 sublevel

    if non_rad_transition_is_collective == 1
      gamma_collective = gamma_zpl + gamma_sideband + gamma_dark;  % total collective decay rate from upper state m=0
    else
      gamma_collective = gamma_zpl + gamma_sideband;  % total collective decay rate from upper state m=0
    end

  else
    gamma_dark = 1/lifetime_dark_1;
    gamma = gamma_zpl + gamma_sideband + gamma_dark;  % total decay rate from upper state m=+1/-1
    gamma_dephasing = 1 / dephasing_lifetime_1;   % decay rate due to local dephasing for m=+1/-1 sublevel

    if non_rad_transition_is_collective == 1
      gamma_collective = gamma_zpl + gamma_sideband + gamma_dark;  % total collective decay rate from upper state m=+1/-1
    else
      gamma_collective = gamma_zpl + gamma_sideband;  % total collective decay rate from upper state m=+1/-1
    end
  end



  last_printed_progress = 0;
  for domain_index = 1:num_domains

    fprintf('\n\nDoing domain %i\n', domain_index);

    Nc = initial_collective_populations(domain_index);   % number of centres initially in collective subspace m=0
    Nnc = initial_noncollective_populations(domain_index);      % number of centres initially in non-collective subspace m=0

    % If there's only one collective atom in the domain, it's basically
    % not collective. While technically the pure collective code will still work
    % the idea of dephasing becomes tricky (there's nowhere for it to go)
    % and dark decays aren't handled well.
    % So in this case, move that one atom to the non-collective population
    if Nc == 1
      Nc = 0;
      Nnc = Nnc + 1;
    end

    fprintf('Collective population    %i\n', Nc);
    fprintf('Noncollective population %i\n', Nnc);

    % We have to zero all the per-domain quantities, since we're going over
    % this loop multiple times

    domain_visible_scattering_rate = zeros(1, num_simulation_time_points);
  
    domain_collective_scattering_rate = zeros(1, num_simulation_time_points);
    domain_visible_collective_scattering_rate = zeros(1, num_simulation_time_points);
  
    domain_noncollective_scattering_rate = zeros(1, num_simulation_time_points);
    domain_visible_noncollective_scattering_rate = zeros(1, num_simulation_time_points);
  
    domain_population_nc = zeros(1, num_simulation_time_points);


    % We have issues if there are no collective atoms at all, i.e. 
    % if Nc = 0, so handle that case separately.

    if Nc == 0
      for time_index = 1 : num_simulation_time_points;
        t = simulation_times(time_index);

        domain_population_nc(time_index) = exp(-gamma*t) * Nnc;
        domain_visible_collective_scattering_rate(time_index) = 0;
        %domain_visible_noncollective_scattering_rate(time_index) = (gamma - gamma_zpl) * domain_population_nc(time_index);
        domain_visible_noncollective_scattering_rate(time_index) = (gamma_sideband) * domain_population_nc(time_index);
        domain_visible_scattering_rate(time_index) = domain_visible_noncollective_scattering_rate(time_index);
      end

      visible_collective_scattering_rate = visible_collective_scattering_rate + domain_visible_collective_scattering_rate;
      visible_scattering_rate = visible_scattering_rate + domain_visible_scattering_rate;
      visible_noncollective_scattering_rate = visible_noncollective_scattering_rate + domain_visible_noncollective_scattering_rate;


%fprintf('location a: \n');
%sum(visible_noncollective_scattering_rate)

      continue; % go the next domain
    end         % end case where a domain has no collective centres

%fprintf('location b')
    % Jmax is the maximum spin subspace. I.e. the one we start in, with Nc0 
    % collective spins. 
    Jmax = Nc/2;

    % The basis vector parameterization we're using is
    % [ Jmax, Jmax  ;  Jmax, Jmax-1  ;  Jmax, Jmax-2  ; ... ;  Jmax, -Jmax  ;
    %    Jmax-1, Jmax-1  ;  Jmax-1, Jmax-2  ;  ...  ;  Jmax-1, -Jmax+1  ;
    %    ...
    %    1/2, 1/2  ;  1/2, -1/2 ]
    %
    % There are (N^2 + 3N) / 2 = Jmax (2Jmax+3) entries in the vector, where N = 2 Jmax
    % Given a specific |J, M>, the corresponding index address in this vector
    % is given by 
    %
    % index = (Jmax - J) (2Jmax + 2J + 3) + J - M + 1 
    %
    % Choosing the highest index which has J=1/2, M=-1/2 which
    % gives index = Jmax (2Jmax+3). Good

    basis_size = (Nc^2 + 3*Nc) / 2;

    matrix = zeros(basis_size, basis_size);
    rho = zeros(basis_size, num_simulation_time_points);

    if offresonant == 1
      rho(1,1) = 1;

    elseif equal_populations == 1
      % We assume due to detunings and dephasing we've been driving
      % to an initial state with equal populations in all M-levels

      for Mindex = 1 : 2*Jmax+1
        rho(Mindex, 1) = 1/(2*Jmax+1);
      end

    elseif gaussian_populations == 1
      % gaussian distribution of populations in the M-levels
      for M = -Jmax : Jmax
        Mindex = Jmax + M + 1;
        s = gaussian_population_width * (2*Jmax+1); % one std deviation
        rho(Mindex, 1) = exp(-M^2/(2*s^2));
      end

      % Normalize the density matrix correctly so it has trace of one
      totalMpopulation = sum(rho(:,1));
      rho(:,1) = rho(:,1) / totalMpopulation;      
    else
      % This is the case where we coherently excite the system

      % Calculate the initial state, assuming we start with all spins down
      % (i.e. in the ground state) and then apply an on-resonant coupling
      % for a specific amount of time, and that we stay only in the J = Jmax
      % manifold
      % That is, H = Omega Jx
      % Take Omega = 1 and choose time such that we get the inversion we want

      % Note that initial_state(1) is M = +Jmax
      % and initial_state(Nc+1) is M = -Jmax
      
      initial_state = zeros(Nc+1, 1);
      initial_state(Nc+1) = 1;

      H = zeros(2*Jmax+1, 2*Jmax+1);

      for M = Jmax : -1 : -Jmax+1
        row = Jmax-M+2;
        column = Jmax-M+1;
        H(row, column) = 0.5*sqrt(Jmax*(Jmax+1) - M*(M-1));
      end

      for M = Jmax-1 : -1 : -Jmax
        row = Jmax-M;
        column = Jmax-M+1;
        H(row, column) = 0.5*sqrt(Jmax*(Jmax+1) - M*(M+1));
      end

      % with 10 atoms takes 3.2 seconds to get full inversion
      % Same with 80 atoms, time to full inversion is independent of atom number
      % Time required for full inversion is t=pi
      % Population in upper state is Nc*sin^2(t/2)
      final_state = zeros(2*Jmax+1, 1);

      time_required = 2*asin(sqrt(proportion_excited_atoms));
      final_state = expm(-i*H*time_required) * initial_state;

      % Take the populations after this coherent driving as the input for
      % the diagonal elements of the density matrix at t=0
      for index = 1 : 2*Jmax+1
        rho(index,1) = (abs(final_state(index))).^2;
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    % end initial condition setup

  
    % matrix is the matrix that describes the Liouvillian
    % for collective decay plus the local atom phase flip process
    % that moves population out of the the |J> collective manifold
    % and into the |J-1/2> collective manifold manifold 
    %
    % d/dt rho = matrix * rho
    %
    % which has the solution rho(t) = exp(matrix * t) * rho(0)
    %

    % This equation is
    %
    % d/dt <J,M| rho | J,M> = g [ (J(J+1) - M(M+1)) rho_J,M+1 - (J(J+1) - M(M-1)) rho_J,M ]
    %                             - 2 J gamma (1 - |M/J|^2) rho_J,M
    %                             + 2(J+1/2) gamma (1 - |(M+1/2) / (J+1/2)|^2) rho_J+1/2,M+1/2
    %                             - (J+M) gamma_isc rho_J,M
    %				  + (J+M+1) gamma_isc rho_J+1/2,M+1/2
    %
    % where the first line is the collective decay, the second line is
    % the dephasing due to local atom pi phase flips causing loss out 
    % of the |J> subspace into the |J-1/2> subspace and
    % the third is the gain into the |J> subspace from the dephasing
    % out of the |J+1/2> subspace
    % and the fourth is the loss from via the isc rate into the shelving levels;
    % we never see those photons so it's a pure loss (note J+M is the total
    % number of excitations).

    % Construct matrix

    for J = Jmax : -1/2 : 1/2
      for M = J : -1 : -J
        row = (Jmax - J) * (2*Jmax + 2*J + 3) + J - M + 1; 

        if M+1 <= J
          column = (Jmax - J) * (2*Jmax + 2*J + 3) + J - (M+1) + 1;
          matrix(row, column) = matrix(row, column) + gamma_collective * (J*(J+1)-M*(M+1));
        end

        column = (Jmax - J) * (2*Jmax + 2*J + 3) + J - M + 1;
        matrix(row, column) = matrix(row, column) - gamma_collective * (J*(J+1)-M*(M-1));
        matrix(row, column) = matrix(row, column) - gamma_dephasing * 2 * J * (1 - (M/J)^2);


        % And now the loss from the isc dark decays
        matrix(row, column) = matrix(row, column) - gamma_dark * (J+M);

        % Now the gain from the dephasing loss in the |Jp, Mp> = |J+1/2, M+1/2> state,
        % since losing an excited atom reduces goes from the J subspace to 
        % J-1/2, as well as taking M to M-1/2 since it's an excited atom we 
        %lost
        if J<Jmax
          Jp = J + 1/2;
          Mp = M + 1/2;
          column = (Jmax - Jp) * (2*Jmax + 2*Jp + 3) + Jp - Mp + 1;
          matrix(row, column) = matrix(row, column) + gamma_dephasing * 2 * Jp * (1 - (Mp/Jp)^2);

          % And now the gain from the isc dark decays
          matrix(row, column) = matrix(row, column) + gamma_dark * (J+M);

        end

      end    % end looping over M
    end      % end looping over J


    % Now actually run the simulation
    integral = 0;
    for time_index = 1 : num_simulation_time_points;
      t = simulation_times(time_index);

      % Don't want to simulate the collective effects for longer than necessary
      % given that it requires matrix exponentials and is slow
      if t < collective_effect_cutoff_time
        % Display progress
        progress = round(100*t/collective_effect_cutoff_time);

        if mod(progress, 10) == 0 && progress ~= last_printed_progress
          last_printed_progress = progress;
          fprintf('%i%%   ', progress);
        end

        rho(:, time_index) = expm(matrix*t) * rho(:,1);
      end

      % Calculate the collective scattering rate from each of the
      % individual J subspaces. There are 2 Jmax - 1 such subspaces.
      % They have dimension 2Jmax+1, 2Jmax, 2Jmax-1... 2.

      % Loop over all the J subspaces
      for J = 1/2 : 1/2 : Jmax

        % For each subspace, construct the diagonal matrix giving the
        % collective scattering rate
        collective_rate_subspace = 0;
        for n = 0 : 2*J;
          M = J - n;
          fm = sqrt(J*(J+1) - M*(M-1));
          index = (Jmax - J) * (2*Jmax + 2*J + 3) + J - M + 1;
        
          % When calculating the fluorescence, make sure we don't include the 
          % non-radiative part even if that's collective, since we can't see it.
          collective_rate_subspace = collective_rate_subspace + (gamma_zpl + gamma_sideband) * fm^2 * rho(index, time_index);
        end

        domain_collective_scattering_rate(time_index) = domain_collective_scattering_rate(time_index) + collective_rate_subspace;

      end

      % The dephasing arises from local spin flips of excited atoms. 
      % For each spin flip, we can consider the projection onto the
      % space |J> and find that the probability of leaving the space
      % and going into |J-1> is (1-|M/J|^). If this occurs we have
      % lost an |e> atom from the collective emsemble and gained one
      % into the non-collective population.
      %
      % The rate at which this occurs is (see p. 174 of notebook)
      %
      % sum_{J=1/2, 1, 3/2, ... Jmax} sum_{-J<=M<=J} 2 J gamma (1-|M/J|^2)
      %
      % I will need the integral of this loss rate from the collective
      % to non-collective manifolds, weighted by exponential decays,
      % so calculate those here.

      collective_loss_rate = 0;

      % Loop over all the J subspaces
      for J = 1/2 : 1/2 : Jmax
        for n = 0 : 2*J;
          M = J - n;
          index = (Jmax - J) * (2*Jmax + 2*J + 3) + J - M + 1;
          collective_loss_rate = collective_loss_rate + 2*J*gamma_dephasing*(1-(M/J)^2) * rho(index, time_index);
        end
      end

      integral = integral + exp(gamma*t) * collective_loss_rate * dt;

      domain_population_nc(time_index) = exp(-gamma*t) * Nnc + exp(-gamma*t) * integral;

      % The experiment can't see photons from the upper levels to the zero-phonon line of the lower
      % levels due to filtering. This means they only measure the rates from
      % the upper levels going into the phonon sidebands, *not* the zero phonon line
      % But the number of photons going into the zpl should be small, so don't break
      % it out explicitly for now
      % The final fluorescence rate is given by 
      %
      % fluorescence rate = visible_collective_scattering_rate + gamma_sideband * non_collective_population
      %
      % To get the visible collective
      % scattering rate we need to scale by the branching ratios - we only see fluorescence
      % into the phonon sidebands, and not into the zpl or the isc states.

      domain_visible_collective_scattering_rate(time_index) = (gamma_sideband/gamma_collective) * domain_collective_scattering_rate(time_index);
      %domain_visible_noncollective_scattering_rate(time_index) = (gamma - gamma_zpl) * domain_population_nc(time_index);
      domain_visible_noncollective_scattering_rate(time_index) = (gamma_sideband) * domain_population_nc(time_index);

      domain_visible_scattering_rate(time_index) = domain_visible_collective_scattering_rate(time_index) + domain_visible_noncollective_scattering_rate(time_index);

    end % end the simulation looping over time points

    % Have reached the end of the simulation a single domain, so add the scattering
    % for this domain to the running total

    visible_collective_scattering_rate = visible_collective_scattering_rate + domain_visible_collective_scattering_rate;
    visible_scattering_rate = visible_scattering_rate + domain_visible_scattering_rate;
    visible_noncollective_scattering_rate = visible_noncollective_scattering_rate + domain_visible_noncollective_scattering_rate;

    %fprintf('location c\n')
  end % end looping over domains

  %fprintf('location d\n')
  % Have finished doing a magnetic sublevel, so add its fluorescence to the total
  total_visible_scattering_rate = total_visible_scattering_rate + visible_scattering_rate;
  total_visible_collective_scattering_rate = total_visible_collective_scattering_rate + visible_collective_scattering_rate;
  total_visible_noncollective_scattering_rate = total_visible_noncollective_scattering_rate + visible_noncollective_scattering_rate;

  %sum(total_visible_noncollective_scattering_rate)*(simulation_times(2)-simulation_times(1))

%sum(total_visible_noncollective_scattering_rate)

end % end looping over magnetic sublevels


% The APD has a finite response time. Need to convolve my simulation result with the detector response.
apd_response_times = APDResponse_16ps(:, 1) * 1e-9;
apd_response_counts = APDResponse_16ps(:, 2);

%test
%apd_response_counts = 1e7*exp(-apd_response_times.^2 / (2*25e-12^2))

apd_response_times = apd_response_times - apd_response_times(1);
adp_norm = sum(apd_response_counts);
apd_response_counts = apd_response_counts / adp_norm;

num_apd_points = length(apd_response_times);


% I need to evaluate s_new(t) = int r(tau) s(t-tau) d tau
% where r(tau) is response function of APD and s(t) is the fluorescence signal
% my simulation gives. s_new is the signal with the detector response
% convolved into it.

total_visible_scattering_rate_convolved = zeros(1, length(total_visible_scattering_rate));
simulation_timestep = simulation_times(2) - simulation_times(1);
response_timestep = apd_response_times(2) - apd_response_times(1);

for time_index = 1:length(total_visible_scattering_rate)
  total = 0;

  for n = 1:num_apd_points

    tau = apd_response_times(n);

    % Convert tau to number of time points in my simulation
    tau_in_simulation_points = floor(tau / simulation_timestep);
    
    if (time_index - tau_in_simulation_points) > 0
      total = total + apd_response_counts(n) * total_visible_scattering_rate(time_index - tau_in_simulation_points);
    end

  end

  total_visible_scattering_rate_convolved(time_index) = total;

end

% Normalise the simulation fluorescence so that the max rate is the same
% for simulation and experimental data
peak_rate_experiment = max(fluorescence_counts_experiment(diamond_number, :));
peak_rate_simulation = max(total_visible_scattering_rate);
peak_rate_simulation_convolved = max(total_visible_scattering_rate_convolved);
amplitude = peak_rate_experiment / peak_rate_simulation;
amplitude_convolved = peak_rate_experiment / peak_rate_simulation_convolved;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Calculate the fit, assumes simulation sample rate is same as experimental sample rate (16ps)
% Do least squares from the peak to 2345 time points later so can compare with other fits

% Find the position of the peak for the convolved theory curve
peak_location_theory = 0;
peak_height_theory = 0;
for index = 1:length(total_visible_scattering_rate_convolved);
  if total_visible_scattering_rate_convolved(index) > peak_height_theory
    peak_location_theory = index;
    peak_height_theory = total_visible_scattering_rate_convolved(index);
  end
end

% Find the position of the peak for the experimental data
peak_location_exp = 0;
peak_height_exp = 0;
for index = 1:length(fluorescence_counts_experiment(diamond_number, :));
  if fluorescence_counts_experiment(diamond_number, index) > peak_height_exp
    peak_location_exp = index;
    peak_height_exp = fluorescence_counts_experiment(diamond_number, index);
  end
end

fluoro_exp = fluorescence_counts_experiment(diamond_number, :) / peak_height_exp;
fluoro_theory = total_visible_scattering_rate_convolved / peak_height_theory;

num_points_to_compare = 2345;

% Calculate mean squares error of fit from peak to end
fit_error = 0;

for time_index = 0:num_points_to_compare-1
  experiment_value = fluoro_exp(peak_location_exp + time_index);
  simulation_value = fluoro_theory(peak_location_theory + time_index);
  fit_error = fit_error + (experiment_value - simulation_value)^2;
end

fit_error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



% Since some photons aren't seen due to filtering (i.e. they come out on the 
% ZPL) we have to back calculate the real number of photons from the visible number
total_collective_photons = gamma_collective/gamma_sideband * sum(total_visible_collective_scattering_rate(1:num_simulation_time_points))*dt;
total_non_collective_photons = (gamma_sideband+gamma_zpl)/gamma_sideband * sum(total_visible_noncollective_scattering_rate(1:num_simulation_time_points))*dt;
total_photons = total_collective_photons + total_non_collective_photons;


fprintf('total collective photons emitted %f\n', total_collective_photons);
fprintf('total noncollective photons emitted %f\n', total_non_collective_photons);
fprintf('total photons emitted %f\n\n', total_photons);

sim_startpoint = 1;
sim_endpoint = 2500;

exp_startpoint = 1500;
exp_endpoint = 4000;
t = simulation_times(sim_startpoint:sim_endpoint) + 0*origin_shift;

linewidth = 2;

set(0, "defaulttextfontsize", 20)  % title
set(0, "defaultaxesfontsize", 16)  % axes labels

figure(1)
plot(t*1e9, total_visible_collective_scattering_rate(sim_startpoint:sim_endpoint), 'r', 'LineWidth', linewidth)
xlabel("Time (ns)")
ylabel("Fluorescence")
title("N=2, initial pop totally inverted")

%figure(1); print N2_initial_pop_middle_of_Dicke_ladder_m0.pdf

%figure(2)
%plot(t, visible_noncollective_scattering_rate(sim_startpoint:sim_endpoint), 'r', 'LineWidth', linewidth)

%figure(2)
%plot(t - origin_shift, visible_scattering_rate(sim_startpoint:sim_endpoint), 'r', 'LineWidth', linewidth)

% Deliberately break out of the code since we're only plotting the theory curves
asdf

amp = 1;

figure(2)
plot(t_experiment(exp_startpoint:exp_endpoint) - origin_shift, fluorescence_counts_experiment(diamond_number, (exp_startpoint:exp_endpoint))/peak_rate_experiment, 'b', 'LineWidth',linewidth)
hold on
plot(t - origin_shift,  amp*amplitude_convolved*total_visible_scattering_rate_convolved(sim_startpoint:sim_endpoint)/peak_rate_experiment, 'r', 'LineWidth', linewidth)
hold off

figure(3)
plot(t_experiment(exp_startpoint:exp_endpoint) - origin_shift, log(fluorescence_counts_experiment(diamond_number, (exp_startpoint:exp_endpoint))/peak_rate_experiment), 'b', 'LineWidth',linewidth)
hold on
plot(t - origin_shift,  log(amp*amplitude_convolved*total_visible_scattering_rate_convolved(sim_startpoint:sim_endpoint)/peak_rate_experiment), 'r', 'LineWidth', linewidth)
hold off

experimental_fluorescence = fluorescence_counts_experiment(diamond_number, 1500:3999)/peak_rate_experiment;
theory_fluorescence = amp*amplitude_convolved*total_visible_scattering_rate_convolved(1:2500)/peak_rate_experiment;
times = t_experiment(1500:3999);

t_ex = t_experiment(exp_startpoint:exp_endpoint) - origin_shift;
t_th = t - origin_shift;
fl_ex = fluorescence_counts_experiment(diamond_number, (exp_startpoint:exp_endpoint))/peak_rate_experiment;
fl_th_50 = total_visible_scattering_rate_convolved(sim_startpoint:sim_endpoint)/peak_rate_experiment;

save('n50_0.4_10_1.8.mat', 't_ex', 't_th', 'fl_ex', 'fl_th_50');

%save('lifetime_plot_data_NV_01.mat', 'experimental_fluorescence', 'theory_fluorescence', 'times');
%save('lifetime_plot_data_NV_17.mat', 'experimental_fluorescence', 'theory_fluorescence', 'times');
%save('lifetime_plot_data_NV_29.mat', 'experimental_fluorescence', 'theory_fluorescence', 'times');
%save('lifetime_plot_data_NV_56.mat', 'experimental_fluorescence', 'theory_fluorescence', 'times');
%save('lifetime_plot_data_NV_86_exp_data_only.mat', 'experimental_fluorescence', 'times');

% 35 centres took 910 seconds
% 40 centres took 1882 seconds
% 45 centres took 3958 seconds
% 50 centres took 6782 seconds

toc


adsf

% Note standard rates for bulk:
% ZPL + sideband: 83 MHz (12ns)
% ISC 0         : 15 MHz
% ISC 1         : 59 MHz
% ZPL + sideband + ISC0:  98 MHz
% ZPL + sideband + ISC1:  142 MHz

amp = 3280;
g0 = 5e7;
g1 = 60e7;
p0 = 0.5;
p1 = 1-p0;
offset = 26.5e-9

figure(1)
plot(t_experiment-offset, fluorescence_counts_experiment(diamond_number, :), 'b', t_experiment, amp*(p0*exp(-g0*t_experiment) + p1*exp(-g1*t_experiment)), 'r');

figure(2)
plot(t_experiment-offset, log(fluorescence_counts_experiment(diamond_number, :)), 'b', t_experiment, log(amp*(p0*exp(-g0*t_experiment) + p1*exp(-g1*t_experiment))), 'r');



asff

for index = 1:num_datasets
  peak(index) = max(fluorescence_counts_experiment(index, :));
end

startpoint = 1;
endpoint = 3600;


% Find the onset shifts by looking at the time at which the
% fluorescence is halfway to the peak, relative to the 
% time at which this happens for the first dataset
% To get what the offset should be do linear interpolation between
% the first point above 0.5 and the one immediately before

onset_shift = zeros(num_datasets, 1);

for time_index = 1:number_of_time_bins
    if fluorescence_counts_experiment(1, time_index) / peak(1) > 0.5
      t1 = t_experiment(time_index-1);
      t2 = t_experiment(time_index);
      h1 = fluorescence_counts_experiment(1, time_index-1) / peak(1);
      h2 = fluorescence_counts_experiment(1, time_index) / peak(1);
      reference_onset_time = t1 + (t2-t1) * (0.5-h1)/(h2-h1);
      break;
    end
end

for dataset_index = 2: num_datasets
  for time_index = 1:number_of_time_bins
    if fluorescence_counts_experiment(dataset_index, time_index) / peak(dataset_index) > 0.5
      t1 = t_experiment(time_index-1);
      t2 = t_experiment(time_index);
      h1 = fluorescence_counts_experiment(dataset_index, time_index-1) / peak(dataset_index);
      h2 = fluorescence_counts_experiment(dataset_index, time_index) / peak(dataset_index);
      onset_shift(dataset_index) = reference_onset_time - (t1 + (t2-t1) * (0.5-h1)/(h2-h1));
      break;
    end
  end
end


% Powers: % 1, 2, 5, 10, 20, 30, 40, 50uW
% Diamonds: NV01, 13, 14, 17, 36

colour_list = {'red' 'green' 'blue' 'cyan' 'magenta' 'yellow' 'black' 'r--' 'g--'};

%All decays get faster as we increase power, effect is greatest with NV01, NV17

%dataset_list = [1, 2, 3, 4, 5, 6, 7, 8];  % NV01 
%dataset_list = [9, 10, 11, 12, 13, 14, 15, 16];  % NV13
%dataset_list = [17, 18, 19, 20, 21, 22, 23, 24];  % NV14
dataset_list = [25, 26, 27, 28, 29, 30, 31, 32];  % NV17
%dataset_list = [33, 34, 35, 36, 37, 38, 39, 40];  % NV36

close all

figure(10)
hold on
for index = 1:length(dataset_list)
  dataset_id = dataset_list(index);
  plot(t_experiment(startpoint:endpoint) + onset_shift(dataset_id), fluorescence_counts_experiment(dataset_id, startpoint:endpoint)/peak(dataset_id), char(colour_list(index)), 'LineWidth', linewidth)
end

figure(2)
hold on
for index = 1:length(dataset_list)
  dataset_id = dataset_list(index);
  plot(t_experiment(startpoint:endpoint) + onset_shift(dataset_id), log(fluorescence_counts_experiment(dataset_id, startpoint:endpoint)/peak(dataset_id)), char(colour_list(index)), 'LineWidth', linewidth)
end


