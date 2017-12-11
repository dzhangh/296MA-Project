% Mini Project for 296MA

% Eco-Solar: Keyu Wang, Haixu Yin, David Zhang

% 1. Find the optimal tilted angle for the solar panel for each month.
% 2. Find the optimal tilted angle for each year.
% 3. Find the corresponding energy generation and payback period.


%% Input all Parameters for each Month

%  Monthly K, optical depth and C, sky diffusion factor and Cn
% ro and I would be constant for the whole year.
K = [0.142,0.144,0.156,0.180,0.196,0.205,0.207,0.201,0.177,0.160,0.149,0.142];
C = [0.058,0.06,0.071,0.097,0.121,0.134,0.136,0.122,0.092,0.073,0.063,0.057];
Cn = [0.48,0.53,0.53,0.58,0.61,0.62,0.67,0.64,0.65,0.59,0.5,0.46];
ro = 0.2;
I = 1300;

% Set up the incremental tilted angle(deg)
beta = 1:1:90;

% Set up the vector for the maximum total solar radiation for each month
Ic_max = zeros(1,12);

% Set up the vector for the optimal angle for each month
beta_opt = zeros(1,12);

% Finding the optimal angle and the corresponding radiation for each month

%% For January
alt_angle = [0,10,20,30,20,10,0];
azi_angle = [63,53,38.6,0,-38.6,-53,-63];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(1)*I*exp(-K(1)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(1)*Cn(1)*I*exp(-K(1)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(1)*I*exp(-K(1)/sind(alt_angle(j)))*(C(1)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_1 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_1(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(1) = max(Ic_avg_1);
beta_opt(1) = find(Ic_avg_1==Ic_max(1));

%% For Feburary
alt_angle = [0,10,20,30,38.9,30,20,10,0];
azi_angle = [75,66,55,38.8,0,-38.8,-55,-66,-75];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(2)*I*exp(-K(2)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(2)*Cn(2)*I*exp(-K(2)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(2)*I*exp(-K(2)/sind(alt_angle(j)))*(C(2)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_2 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_2(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(2) = max(Ic_avg_2);
beta_opt(2) = find(Ic_avg_2==Ic_max(2));

%% For March
alt_angle = [0,10,20,30,40,50,40,30,20,10,0];
azi_angle = [90,81.6,72.5,60.8,45,0,-45,-60.8,-72.5,-81.6,-90];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(3)*I*exp(-K(3)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(3)*Cn(3)*I*exp(-K(3)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(3)*I*exp(-K(3)/sind(alt_angle(j)))*(C(3)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_3 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_3(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(3) = max(Ic_avg_3);
beta_opt(3) = find(Ic_avg_3==Ic_max(3));

%% For April
alt_angle = [0,10,20,30,40,50,60,60,50,40,30,20,10,0];
azi_angle = [104,95.7,87.4,78,67.5,51.7,15.2,-15.2,-51.7,-67.5,-78,-87.4,-95.7,-104];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(4)*I*exp(-K(4)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(4)*Cn(4)*I*exp(-K(4)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(4)*I*exp(-K(4)/sind(alt_angle(j)))*(C(4)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_4 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_4(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(4) = max(Ic_avg_4);
beta_opt(4) = find(Ic_avg_4==Ic_max(4));

%% For May
alt_angle = [0,10,20,30,40,50,60,70,60,50,40,30,20,10,0];
azi_angle = [116,107,99,91,82.3,71.3,55.3,0,-55.3,-71.3,-82.3,-91,-99,-107,-116];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(5)*I*exp(-K(5)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(5)*Cn(5)*I*exp(-K(5)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(5)*I*exp(-K(5)/sind(alt_angle(j)))*(C(5)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_5 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_5(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(5) = max(Ic_avg_5);
beta_opt(5) = find(Ic_avg_5==Ic_max(5));

%% For June
alt_angle = [0,10,20,30,40,50,60,70,73.3,70,60,50,40,30,20,10,0];
azi_angle = [120,112.5,104.5,96.5,88.5,79,65.6,37.4,0,-37.4,-65.6,-79,-88.5,-96.5,-104.5,-112.5,-120];

Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(6)*I*exp(-K(6)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(6)*Cn(6)*I*exp(-K(6)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(6)*I*exp(-K(6)/sind(alt_angle(j)))*(C(6)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end

Ic_avg_6 = zeros(length(beta),1);

for i = 1:length(beta)
    Ic_avg_6(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(6) = max(Ic_avg_6);
beta_opt(6) = find(Ic_avg_6==Ic_max(6));

%% For July
alt_angle = [0,10,20,30,40,50,60,70,60,50,40,30,20,10,0];
azi_angle = [116,107,99,91,82.3,71.3,55.3,0,-55.3,-71.3,-82.3,-91,-99,-107,-116];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(7)*I*exp(-K(7)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(7)*Cn(7)*I*exp(-K(7)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(7)*I*exp(-K(7)/sind(alt_angle(j)))*(C(7)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_7 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_7(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(7) = max(Ic_avg_7);
beta_opt(7) = find(Ic_avg_7==Ic_max(7));

%% For August
alt_angle = [0,10,20,30,40,50,60,60,50,40,30,20,10,0];
azi_angle = [104,95.7,87.4,78,67.5,51.7,15.2,-15.2,-51.7,-67.5,-78,-87.4,-95.7,-104];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(8)*I*exp(-K(8)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(8)*Cn(8)*I*exp(-K(8)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(8)*I*exp(-K(8)/sind(alt_angle(j)))*(C(8)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_8 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_8(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(8) = max(Ic_avg_8);
beta_opt(8) = find(Ic_avg_8==Ic_max(8));

%% For September
alt_angle = [0,10,20,30,40,50,40,30,20,10,0];
azi_angle = [90,81.6,72.5,60.8,45,0,-45,-60.8,-72.5,-81.6,-90];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(9)*I*exp(-K(9)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(9)*Cn(9)*I*exp(-K(9)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(9)*I*exp(-K(9)/sind(alt_angle(j)))*(C(9)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_9 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_9(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(9) = max(Ic_avg_9);
beta_opt(9) = find(Ic_avg_9==Ic_max(9));

%% For October
alt_angle = [0,10,20,30,38.9,30,20,10,0];
azi_angle = [75,66,55,38.8,0,-38.8,-55,-66,-75];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(10)*I*exp(-K(10)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(10)*Cn(10)*I*exp(-K(10)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(10)*I*exp(-K(10)/sind(alt_angle(j)))*(C(10)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_10 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_10(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(10) = max(Ic_avg_10);
beta_opt(10) = find(Ic_avg_10==Ic_max(10));

%% For November
alt_angle = [0,10,20,30,20,10,0];
azi_angle = [63,53,38.6,0,-38.6,-53,-63];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(11)*I*exp(-K(11)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(11)*Cn(11)*I*exp(-K(11)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(11)*I*exp(-K(11)/sind(alt_angle(j)))*(C(11)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_11 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_11(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(11) = max(Ic_avg_11);
beta_opt(11) = find(Ic_avg_11==Ic_max(11));

%% For December
alt_angle = [0,10,20,26.6,20,10,0];
azi_angle = [58.7,47.7,31.2,0,-31.2,-47.7,-58.7];

% Set up the matrix for each term
Ic_1 = zeros(length(beta),length(alt_angle));
Ic_2 = zeros(length(beta),length(alt_angle));
Ic_3 = zeros(length(beta),length(alt_angle));
Ic = zeros(length(beta),length(alt_angle));

for i = 1:length(beta)
    for j = 1:length(alt_angle)
        Ic_1(i,j) = Cn(12)*I*exp(-K(12)/sind(alt_angle(j)))*cosd(alt_angle(j))*(cosd(azi_angle(j))*sind(beta(i))+sind(alt_angle(j))*cosd(beta(i)));
        Ic_2(i,j) = C(12)*Cn(12)*I*exp(-K(12)/sind(alt_angle(j)))*cosd(beta(i)/2)^2;
        Ic_3(i,j) = ro*Cn(12)*I*exp(-K(12)/sind(alt_angle(j)))*(C(5)+sind(alt_angle(j)))*sind(beta(i)/2)^2;
        Ic(i,j) = Ic_1(i,j) + Ic_2(i,j) + Ic_3(i,j);
        j = j + 1;
    end
    i = i + 1;
end
Ic_avg_12 = zeros(length(beta),1);

% Find the average solar radiation for this month
for i = 1:length(beta)
    Ic_avg_12(i) = mean(Ic(i,:));
end

% Find the maximum radiation and the corresponding tilted angle
Ic_max(12) = max(Ic_avg_12);
beta_opt(12) = find(Ic_avg_12==Ic_max(12));


%% Making the plots for max radiation and optimal angle

month = 1:12;

% plot(month,Ic_max,'g-*')
% xlabel('Month(Jan-Dec)')
% ylabel('Total Radiation(W/m^2)')
% title('Maximum Total Radiation on each Month')
% grid on


% plot(month,beta_opt,'r--o')
% xlabel('Month(Jan-Dec)')
% ylabel('Optimal Angle(Degree)')
% title('Optimal Tilted Angle for the Panels on each Month')
% grid on

%% Finding the optimal angle for the year
Ic_year = zeros(length(beta),1);
for i = 1:length(beta)
    Ic_year(i)=Ic_avg_1(i)+Ic_avg_2(i)+Ic_avg_3(i)+Ic_avg_4(i)+Ic_avg_5(i)+Ic_avg_6(i)+Ic_avg_7(i)+Ic_avg_8(i)+Ic_avg_9(i)+Ic_avg_10(i)+Ic_avg_11(i)+Ic_avg_12(i);
end

Ic_max_year = max(Ic_year);

% The optimal tilted angle for a year
beta_opt_year = find(Ic_year==Ic_max_year);

% Comparing monthly radiation with set tilted angle and changing tilted
% angle

% Max radiation for set tilte angle
Ic_max_month = Ic_max_year/12;
% Max radiation for changing tilted angle
Ic_max_opt = mean(Ic_max);

%Average monthly Ic vs Tilt Angle
Ic_avg = Ic_year/12;

% plot(beta,Ic_avg,'g-*')
% xlabel('Tilt Angle [deg]')
% ylabel('Average Monthly Solar Radiation(W/m^2)')
% title('Average Monthly Solar Radiation vs. Tilt Angle')
% grid on

%% Find the Energy Generated per month

% Input the hours of the day for each month
hrs = [10,11,12,13.5,14.5,14.75,14.5,13.5,12.25,11,10,9.5];

% Input the number of days of each month
num_days = [31,28,31,30,31,30,31,31,30,31,30,31];

% Find the radiation for a tilted angle of 45 degree for each month;
Ic_45 = [Ic_avg_1(45),Ic_avg_2(45),Ic_avg_3(45),Ic_avg_4(45),Ic_avg_5(45),Ic_avg_6(45),Ic_avg_7(45),Ic_avg_8(45),Ic_avg_9(45),Ic_avg_10(45),Ic_avg_11(45),Ic_avg_12(45)];

% plot(month,Ic_45,'g-*')
% xlabel('Month (Jan-Dec.)')
% ylabel('Solar Radiation(W/m^2)')
% title('Solar Radiation at the Optimal PV Tilt Angle (45 deg)')
% grid on

% With a panel of area 1m^2 and efficiency being 16%, calculate the energy generated for each month
A = 1; %m^2
W_max = 100; %Rated maximum output wattage
eff = 0.16; %Conventional PV
W_month = zeros(length(hrs),1);
E_month = zeros(length(hrs),1);
for i = 1:length(hrs)
    W_month(i) = Ic_45(i)*A*eff;
    if W_month(i)> W_max
        E_month(i) = W_max*hrs(i)*num_days(i)/1000;
    else
        E_month(i) = W_month(i)*hrs(i)*num_days(i)/1000; %kWh/month
    end
end

% plot(month,W_month,'g-*')
% xlabel('Month (Jan-Dec)')
% ylabel('Power Output (W)')
% title('Power Output at the Optimal PV Tilt Angle (45 deg)')
% grid on

% Energy generated for each year
E_year = sum(E_month); %kWh/year

% The amount of utility saved
earn_year = 0.16*E_year; %assume 16 cents / kWh

% Calculate the paypack period
cost = 350; %Cost of Eco-Solar in $
pay_period = cost/earn_year; %Estimated payback period

