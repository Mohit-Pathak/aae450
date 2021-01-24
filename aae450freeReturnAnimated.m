% AAE 532 Orbit Mechanics
% Purdue University
% HW 9, Problem 1
% Fall 2020 - 11/13/2020
% By: Mohit Pathak

clearvars, close all
clc

%%
rMoon = 1738.2;
muMoon = 4902.8005821478;
aMoon = 384400;
rEarth = 6378.1363;
muEarth = 398600.4415;

%% pt. a
alpha1 = 38.2529;
gamPlus1 = 33.7046;
eta = 180-alpha1;
beta2 = 180-eta-gamPlus1;
alpha2 = 180-beta2;

%% pt. b
TA = 173.8;
altEarth = 190;
rp = rEarth + altEarth;
eT = (rp-aMoon)/(aMoon*cosd(TA)-rp);
aT = rp/(1-eT);
ra = aT*(1+eT);
IP = 2*pi*sqrt(aT^3/muEarth);
specEnergy = -muEarth/(2*aT);
E = rad2deg(acos((aT-aMoon)/(aT*eT)));
M = deg2rad(E) - eT*sind(E);
n = (2*pi)/IP;
tMinustP = M/n;
n2 = sqrt(muEarth/(aMoon^3));
phi = rad2deg(deg2rad(TA) - (n2*tMinustP));

%% pt. c
vMinus = sqrt((2*muEarth/aMoon)-(muEarth/aT));
vPlus = vMinus;
rPlus = aMoon;
p = aT*(1-eT^2);
h = sqrt(muEarth*p);
gamMinus = rad2deg(acos(h/(rPlus*vPlus)));
vMoon = sqrt(muEarth/aMoon);
vInfMoon = sqrt((vMoon^2)+(vMinus^2)-(2*vMoon*vMinus*cosd(gamMinus)));
flybyAng = rad2deg(2*acos(((vMinus^2)-(vMoon^2)-(vInfMoon^2))/(-2*vMoon*vInfMoon)));
eH = 1/sind(flybyAng/2);
specEnergyH = 0.5*(vInfMoon^2);
aH = -muMoon/(vInfMoon^2);
rpPass = aH*(1-eH);
passAlt = rpPass-rMoon;
dvEq = sqrt((vInfMoon^2)+(vInfMoon^2)-(2*vInfMoon*vInfMoon*cosd(flybyAng)));
beta1 = (180-2*gamMinus)/2;
alpha1 = 180-beta1; % alpha negative because s/c appraoching Earth

%% pt. d plotting
% pt. i) earth centered frame
% plot earth orbit
% subplot(2,2,[3,4])
% suptitle('Lunar Free-Return Trajectories');
% figure()
suptitle('Lunar Free Return Trajectories - By: Mohit Pathak');
subplot(2,1,1);
thetaStar = (linspace(0,2*pi, 3600000))';
xEarth = rEarth.*cos(thetaStar);
yEarth = rEarth.*sin(thetaStar);
plot(xEarth,yEarth, 'g', 'LineWidth', 2)
hold on;
grid on
% plot moon orbit
thetaStar = (linspace(0,2*pi, 3600000))';
rTemp2 = aMoon;
xMoon = rTemp2.*cos(thetaStar);
yMoon = rTemp2.*sin(thetaStar);
hold on;
plot(xMoon,yMoon, 'k--', 'LineWidth', 1)
axis equal
hold on;
plot(xMoon(1738000),yMoon(1738000),'k*','LineWidth', 2)
grid on
% plot parking orbit around Earth
rTemp3 = rp;
xParking = rTemp3.*cos(thetaStar);
yParking = rTemp3.*sin(thetaStar);
hold on;
plot(xParking,yParking, 'b', 'LineWidth', 2)
hold on;
grid on
axis equal
% plot transfer arcs
% plot out-bound arc
thetaStarOutbound = (linspace(0,(173.8*pi/180)))';
domega = deg2rad(-12.4);
rTemp3 = p./(1+eT.*cos(thetaStarOutbound));
xtransferOut = rTemp3.*cos(thetaStarOutbound);
ytransferOut = rTemp3.*sin(thetaStarOutbound);
hold on;

% subplot(2,2,1)
% plot(xEarth,yEarth, 'g', 'LineWidth', 2)

% animate orbit
curve = animatedline('LineWidth',2);
% plot(xtransferOut,ytransferOut, 'b-', 'LineWidth', 1)
for i=1:length(xtransferOut)
    hold on
    addpoints(curve,xtransferOut(i), ytransferOut(i));
    head = scatter(xtransferOut(i), ytransferOut(i),'filled','MarkerFaceColor','b');
    drawnow
    % pause(0.01);
    delete(head);
    if(i == 93)
       hold on;
       plot(xtransferOut(93),ytransferOut(93), 'b-<', 'LineWidth', 1)
    end
%     subplot(2,2,1);
%     hold on;
%     xlim([-90000,90000]); ylim([-60000, 60000]);
%     addpoints(curve,xtransferOut(i), ytransferOut(i));
%     head = scatter(xtransferOut(i), ytransferOut(i),'filled','MarkerFaceColor','b');
%     drawnow
%     pause(0.01);
%     delete(head);
end
hold on
grid on
axis equal
% plot in-bound arc
thetaStarInbound = (linspace(-(173.8*pi/180),0))';
rTemp3 = p./(1+eT.*cos(thetaStarInbound));
xtransferOut = rTemp3.*cos(thetaStarInbound + domega);
ytransferOut = rTemp3.*sin(thetaStarInbound + domega);
hold on;
curve2 = animatedline('LineWidth',2, 'Color', 'g');
% plot(xtransferOut,ytransferOut, 'g-', 'LineWidth', 1)
hold on;
for i=1:length(xtransferOut)
    addpoints(curve2,xtransferOut(i), ytransferOut(i));
    head = scatter(xtransferOut(i), ytransferOut(i),'filled','MarkerFaceColor','g');
    drawnow
    % pause(0.01);
    delete(head);
    if(i == 7)
       hold on;
       plot(xtransferOut(7),ytransferOut(7), 'g->', 'LineWidth', 1)
    end
end
hold on;
grid on
% make everything look nice :)
title('Full Trajectory View');
xlabel('Distance x (km)');
ylabel('Distance y (km)');
axis equal
% plot moon
% figure()
subplot(2,1,2)
axis equal
xlim([-10000,7000]); ylim([-7000, 7000]);
title("Moon-Centered View")
hold on
thetaStar = (linspace(0,2*pi))';
rTemp3 = rMoon;
xtransferOut = rTemp3.*cos(thetaStar);
ytransferOut = rTemp3.*sin(thetaStar);
hold on;
plot(xtransferOut,ytransferOut, 'k-', 'LineWidth', 2);
fill(xtransferOut,ytransferOut, 'k-')
hold on;
hold on;
grid on
axis equal
% plot hyperbolic orbit
thetaStar = (linspace(-pi,pi))';
pH = aH*(1-eH^2);
rTemp4 = pH./(1+eH.*cos(thetaStar));
xtransferOut = rTemp4.*cosh(thetaStar);
ytransferOut = rTemp4.*sinh(thetaStar);
xtransferOut = xtransferOut-((2*rpPass));
curve3 = animatedline('LineWidth',2);
hold on;
for i=1:length(xtransferOut)
    addpoints(curve3,xtransferOut(i), ytransferOut(i));
    head = scatter(xtransferOut(i), ytransferOut(i),'filled','MarkerFaceColor','r');
    drawnow
    pause(0.01);
    delete(head);
end
%plot(xtransferOut-((2*rpPass)),ytransferOut, 'r-', 'LineWidth', 2);
hold on;
grid on


