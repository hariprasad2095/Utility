clear all
close all
prompt='Base station 1 co-ordinates';
BS1=[0 0 0];
bs1_X = BS1(1);
bs1_Y = BS1(2);
bs1_Z = BS1(3);
prompt='Base station 2 co-ordinates';
BS2=[10,10,0];
bs2_X = BS2(1);
bs2_Y = BS2(2);
bs2_Z = BS2(3);
prompt='UE_co-ordinates co-ordinates';
UE=[1,1,0];
Ue_X = UE(1);
Ue_Y = UE(2);
Ue_Z = UE(3);
%% Assumed UE speed = 65mph.
UE_speed = 65;
distance_bs1_bs2 = sqrt( (bs1_X-bs2_X)^2 + (bs1_Y-bs2_Y)^2); % in miles
UE_time_taken = distance_bs1_bs2/UE_speed; % in hour
time_interval = 0:1:(UE_time_taken*3.6e6); % in milliseconds
elapsed_time = 0;
UE_position=UE; 
unit_step_dist= 1.80556e-5; % miles per millisecond
range_of_bs1 = 0.7*distance_bs1_bs2;
range_of_bs2 = 0.3*distance_bs1_bs2;
UE_pos_arr = [Ue_X Ue_Y];
X1 = bs1_X;
Y1 = bs1_Y;
X2= bs2_X;
Y2 = bs2_Y;
%% Making the UE move
while elapsed_time <= UE_time_taken
    
    %% make the UE move
    % find points that are equidistant. 
    syms x_ue y_ue
    %% [sol_x,sol_y] = solve(x_ue^2 + y_ue^2 -2*X1*x_ue - 2*Y1*y_ue == unit_step_dist - X1^2 - Y1^2 , ...
    %%                       x_ue^2 + y_ue^2 -2*X2*x_ue - 2*Y2*y_ue == distance_bs1_bs2- unit_step_dist - X2^2 - Y2^2 );

    [sol_x,sol_y] = solve(x_ue^2 + y_ue^2 -2*X1*x_ue - 2*Y1*y_ue == unit_step_dist - X1^2 - Y1^2 , ...
                            (x_ue-bs1_X)*(bs2_Y-bs1_Y)-(bs2_X-bs1_X)*(y_ue-bs1_Y)==0);
    % Update the time
    X= double(real(sol_x));
    Y = double(real(sol_y));
    UE_new_pos = [X(1) Y(1)];
    X1=X(1);
    Y1=Y(1);
    UE_pos_arr = [UE_pos_arr UE_new_pos];
    distance_bs1_bs2 = distance_bs1_bs2- unit_step_dist;
    elapsed_time = elapsed_time+0.001;
end




