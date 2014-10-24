%For running the GPR on simulated data, then performing validations and
%evaluating performance of the generated model.
%By: Subhash Gubba
%8/12/2014 - OCCaM Lab - #SHUML

clear all; close all;
addpath('../../..');
[sweep,sweep_vecs] = generate_sweep;                            %Create state space sweep
testRegression_nd_vel;                                              %Run GPR
save('m_vel.mat','mF_vel');
testRegression_nd_accel;
load('m_vel.mat','mF_vel');
hold off;
theoretical_ds = one_step_sim(E,sweep);                         %Generate outputs from sweep
experimental_ds = [mF_vel mF_accel];
error_calc_gpr(sweep,sweep_vecs,theoretical_ds,experimental_ds);