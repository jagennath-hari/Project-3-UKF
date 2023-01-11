clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.01*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;
%% Calculate Kalmann Filter
for i = 1:length(sampledTime)
    %% FILL IN THE FOR LOOP
    acc = sampledData(i).acc; %loading accelearation
    angVel = sampledData(i).omg; %loading angular velocity
    dt = double(sampledTime(i) - prevTime); %Computing dt

    % Prediction Step
    [covarEst, uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt);

    % Measurement Model
    z_t = [transpose(vel(i,:)); transpose(angVel2(i,:))];
    
    % Update Step

    [uCurr, covar_curr] = upd_step(z_t, covarEst, uEst);

    % Changing Prev Time

    prevTime = sampledTime(i);

    % Adding uCurr to Saved State variable
    
    savedStates(:,i) = uCurr;

    % Changing uPrev to uCurr and same for the covariance
    
    uPrev = uCurr;
    covarPrev = covar_curr;
end

plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);