function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    %%

 C = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]; %C matrix
 R = diag([0.0000001 0.0000001 0.0000001 0.0000001 0.0000001 0.0000001]); %R matrix found by tuning

 Kt = (covarEst * transpose(C))*pinv((((C * covarEst * transpose(C)) + R))); %Kalman gain 

uCurr = double(uEst + (Kt * (z_t - (C * uEst)))); %Computing current mean
covar_curr = double(covarEst - (Kt * C * covarEst)); %Computing current covariance

end

