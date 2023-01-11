function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
    %%
%Noise values
ngx = 1;
ngy = 1;
ngz = 1;
nax = 1;
nay = 1;
naz = 1;
bgx = 1;
bgy = 1;
bgz = 1;
bax = 1;
bay = 1;
baz = 1; %Noise values found by tuning

%%
%Noise Matrices
ng = [ngx; ngy; ngz]; %ng matrix
na = [nax; nay; naz]; %na matrix
nbg = [bgx; bgy; bgz]; %nbg matrix
nba = [bax; bay; baz]; %nba matrix


%%
%Compute sigma points
X = []; %Matrix for holding sigma points
n = 15 + 12; %Numer of sigma points 15 states + 12 for augmentation
a = 0.001; %alpha value
k = 1; %kappa value
b = 2; %beta value
lambda = ((a^2)*(n + k)) - n; %lambda value
muAug = [uPrev; zeros([12, 1])]; %Augmented mean
pAug = [covarPrev zeros([15, 12]); zeros([12, 15]) diag([ng; na; nbg; nba])]; %Augmented covariance matrix
covChol =  chol(pAug, "lower"); %Using cholesky decomposition on augmented covariance matrix
for i = 1 : 55 %for loop for all sigma points 55 total points
    if(i == 1) %if statement for computing the first sigma point
        X(1 : 27, i) = muAug; %first sigma point
    elseif(i ~= 1 && i <= 28) %else if statement to compute the positive sigma points
        X(1 : 27, i) = muAug + (sqrt(n + lambda) * covChol(:, i - 1)); %Computing positive sigma points
    elseif(i > 28) %else if statement to compute to negative sigma points
        X(1 : 27, i) = muAug - (sqrt(n + lambda) * covChol(:, rem(i, 28))); %Computing negative sigma points
    end %ending if statement
end %end for loop
%%
%Propogate Sigma Points
Xt = []; %empty matirx holding propogation of sigma points
wm = [angVel(1,1); angVel(2,1); angVel(3,1)]; %Angular velocity
am = [acc(1,1); acc(2,1); acc(3,1)]; %acceleration

for i = 1 : 55 %for loop for all sigma points
    R = [cos(X(5, i))*cos(X(6, i)), cos(X(6, i))*sin(X(4, i))*sin(X(5, i)) - cos(X(4, i))*sin(X(6, i)), sin(X(4, i))*sin(X(6, i)) + cos(X(4, i))*cos(X(6, i))*sin(X(5, i));
    cos(X(5, i))*sin(X(6, i)), cos(X(4, i))*cos(X(6, i)) + sin(X(4, i))*sin(X(5, i))*sin(X(6, i)), cos(X(4, i))*sin(X(5, i))*sin(X(6, i)) - cos(X(6, i))*sin(X(4, i));
    -sin(X(5, i)), cos(X(5, i))*sin(X(4, i)), cos(X(4, i))*cos(X(5, i))]; %Rotation Matrix

    G = pinv([cos(X(5, i))*cos(X(6, i)), -sin(X(6, i)), 0;
    cos(X(5, i))*sin(X(6, i)), cos(X(6, i)), 0;
    -sin(X(5, i)), 0, 1])* R; %Euler Rates

    Xt(1 : 15, i) = X(1 : 15, i) + (dt*[X(7 : 9, i); G * (wm - X(10 : 12, i) - X(16 : 18, i)); [0 ; 0; -9.81] + (R * (am - X(13 : 15, i) - X(19 : 21, i))); X(22 : 24, i); X(25 : 27, i)]); %Propogation of sigma points
end %ending for loop
%%

%Compute Mean
ut =[]; %empty mean matrix
for i = 1 : 55 %for loop for propogated sigma points
    if (i == 1) %if statement for the first sigma point
        Wm = lambda / (n+lambda); %Computing Wm for the first sigma point
        ut = Wm * Xt(:, i); %Computing the mean for the first sigma point
    else %else statement
        Wm = 1 / (2 * (n + lambda)); %Computing Wm for rest of the sigma points
        ut = ut + (Wm * Xt(:, i)); %Computing the mean for rest of the sigma points
    end %end statement
end %end statement
uEst = ut; %returning Estimated mean

%%
%Compute Covariance

sigmt = []; %empty Covariance matrix
for i = 1 : 55 %for loop for propgated sigma points
    if (i == 1) %if statement for the first sigma point
        Wc = (lambda / (n + lambda)) + (1 - a^2 + b); %Computing Wc for the first sigma point
        sigmt = Wc * (Xt(:, i) - uEst) * transpose((Xt(:, i) - uEst)); %Computing the Covariance for the first sigma point
    else %else statment
        Wc = 1/(2 * (n + lambda)); %Computing Wc for rest of the sigma points
        sigmt = sigmt + (Wc * (Xt(:, i) - uEst) * transpose((Xt(:, i) - uEst))); %Computing the Covariance for rest of the sigma points
    end %end statment for if statment
end %end statment for loop
covarEst = sigmt; %returning covariance

end




