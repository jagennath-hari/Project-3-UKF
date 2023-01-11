function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    %%
%Noise values
n1 = 0.04; 
n2 = 0.04;
n3 = 0.04;
n4 = 0.1;
n5 = 0.1;
n6 = 0.1; %Noise values found by tuning
Rt = diag([n1; n2; n3; n4; n5; n6]); %diagonal matrix of noises
covChol = chol(covarEst, "lower"); %cholesky decomposition
a = 0.001; %alpha value
b = 2; %beta value
n = 15; %number of sigma points
k = 1; %kappa value
lambda = ((a^2)*(n + k)) - n; %lambda value
%%
R_cb = [0.707 -0.707 0; -0.707 -0.707 0; 0 0 -1]; %Rotation of camera wrt body frame
t_cb = [-0.04; 0.0; -0.03];  % Translation of camera wrt body frame
T_cb =   [R_cb t_cb;0 0 0 1];  % Transformation matrix of body frame wrt camera frame
T_bc = (T_cb)^-1; % Transformation matrix of camera wrt body
R_bc = T_bc(1:3,1:3); %Rotaion matrix of body in camera frame
adj = [R_bc -R_bc*skew(T_bc(1:3,4));zeros(3,3) zeros(3,3)]; %adjoint matrix
%%
%Compute the sigma points

X = []; %Matrix for holding sigma points
for i = 1 : 31 %for loop for all sigma points
    if(i == 1) %if statement for computing the first sigma point
        X(1 : 15, i) = uEst;  %first sigma point
    elseif(i ~= 1 && i <= 16) %else if statement to compute the positive sigma points
        X(1 : 15, i) = uEst + (sqrt(n + lambda) * covChol(:, i - 1)); %Computing positive sigma points
    elseif(i > 16) %else if statement to compute to negative sigma points
        X(1 : 15, i) = uEst - (sqrt(n + lambda) * covChol(:, rem(i, 16))); %Computing negative sigma points
    end %end statement for if statments
end %end statment for loop

%%
%Propogate the sigma points
Zt = []; %empty matirx holding propogation of sigma points

for i = 1 : 31 %for loop for all sigma points
    R = ([cos(X(5, 1))*cos(X(6, 1)), cos(X(6, 1))*sin(X(4, 1))*sin(X(5, 1)) - cos(X(4, 1))*sin(X(6, 1)), sin(X(4, 1))*sin(X(6, 1)) + cos(X(4, 1))*cos(X(6, 1))*sin(X(5, 1));
    cos(X(5, 1))*sin(X(6, 1)), cos(X(4, 1))*cos(X(6, 1)) + sin(X(4, 1))*sin(X(5, 1))*sin(X(6, 1)), cos(X(4, 1))*sin(X(5, 1))*sin(X(6, 1)) - cos(X(6, 1))*sin(X(4, 1));
    -sin(X(5, 1)), cos(X(5, 1))*sin(X(4, 1)), cos(X(4, 1))*cos(X(5, 1))]); %Rotation matrix
    Zt(1 : 3, i) = (R_cb * transpose(R) * X(7 : 9, i)) - (R_cb * skew(T_cb(1 : 3, 4)) * transpose(R_cb) * z_t(4 : 6, 1)); %Computing linear velocity
end %end statement

%%
%Computing the  zut
zut = []; %empty matrix for holding zut
for i = 1 : 31 %for loop for all sigma points
    if (i == 1) %if statement for the first sigma point
        Wm = lambda / (n+lambda); %Computing Wm for the first sigma point
        zut = Wm * Zt(:, i); %Computing the zut for the first sigma point
    else %else statment
        Wm = 1 / (2 * (n + lambda)); %Computing Wm for rest of the sigma points
        zut = zut + (Wm * Zt(:, i)); %Computing the zut for rest of the sigma points
    end %end for if statement
end %end statment for loop

%%
%Computing the Ct
Ct = []; %empty matrix for Ct
for i = 1 : 31 %for loop for all sigma points
    if (i == 1) %if statement for first sigma point
        Wc = (lambda / (n + lambda)) + (1 - a^2 + b); %Computing Wc for the first sigma point
        Ct = Wc * (X(1 : 15, i) - uEst) * (transpose(Zt(:, i) - zut)); %Computing the Ct for the first sigma point
    else %else statement
        Wc = 1/(2 * (n + lambda)); %Computing Wc for rest of the sigma points
        Ct = Ct + (Wc * (X(1 : 15, i) - uEst) * (transpose(Zt(:, i) - zut))); %Computing the Ct for rest of the sigma points
    end %end statment for if statement
end %end statment for loop
%%
%Computing the St
St = []; %empty matrix for St
for i = 1 : 31 %for loop for all sigma points
    if (i == 1) %if statement for first sigma point
        Wc = (lambda / (n + lambda)) + (1 - a^2 + b); %Computing Wc for the first sigma point
        St = (Wc * (Zt(:, i) - zut) * (transpose(Zt(:, i) - zut))); %Computing the St for the first sigma point 
    else %else statment
        Wc = 1/(2 * (n + lambda)); %Computing Wc for rest of the sigma points
        St = (St + (Wc * (Zt(:, i) - zut) * (transpose(Zt(:, i) - zut)))); %Computing the St for rest of the sigma points
    end %end statement for if condition
end %end statment for loop
St = St + Rt(1:3, 1:3); %Recomputing St by adding noise
%%
%Computing the Filter Gain
Kt = Ct / St; %Filter Gain
%%
%Computing the Filtered Mean and Covariance
uCurr = uEst + (Kt*(z_t(1 : 3, 1) - zut)); %current mean
covar_curr = covarEst - (Kt * St * transpose(Kt)); %current covaraince

end

