%

addpath('../../..');
testRegression_nd;                                              %Run GPR
hold off;
[sweep,sweep_vecs] = generate_sweep;                            %Create state space sweep
theoretical_ds = one_step_sim(E,sweep);                         %Generate outputs from sweep
experimental_ds = 
error_calc_gpr(sweep,sweep_vecs,theoretical_ds,experimental_ds);