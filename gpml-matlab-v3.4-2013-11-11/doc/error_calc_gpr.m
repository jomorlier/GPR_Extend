function err = error_calc_gpr(sweep,sweep_vecs,theoretical_ds,experimental_ds)
%Computes error between the theoretical and experimental results, should be
%made modular with dimensions of input
    
%     Error is the closeness of the change in one timestep
    err = (((theoretical_ds-sweep(:,1:2)) - (experimental_ds-sweep(:,1:2)))./(theoretical_ds-sweep(:,1:2))); %percent error calculation
%     err = (((theoretical_ds-sweep(:,1:2)) - (experimental_ds-sweep(:,1:2))));
    %disp(err);
    err_condition = cond(err); %finds the condition of the matrix
    disp(['condition = ', num2str(err_condition)]);
    disp(['mean-abs-error = ', num2str(mean(abs(err)))]);%the mean absolute error
    
    S = [length(sweep_vecs.theta_vec) length(sweep_vecs.omega_vec) length(sweep_vecs.cntrl_vec)];
    err_mat_p = reshape(err(:,1),S(1),S(2),S(3)); %change position error into 3D matrix
    err_mat_v = reshape(err(:,2),S(1),S(2),S(3)); %change velocity error into 3D matrix
    err_u_p = mean(abs(err_mat_p),3); %find the average along the control signal (u) dimension
    err_u_v = mean(abs(err_mat_v),3);
    
    figure;
    pcolor(sweep_vecs.theta_vec',sweep_vecs.omega_vec', err_u_p');
    title('Position vs. Velocity vs. Error in Change in Position');
    xlabel('Position');
    ylabel('Velocity');
    colorbar;
    figure;
    pcolor(sweep_vecs.theta_vec',sweep_vecs.omega_vec', err_u_v');
    title('Position vs. Velocity vs. Error in Change in Velocity');
    xlabel('Position');
    ylabel('Velocity');
    colorbar;
end