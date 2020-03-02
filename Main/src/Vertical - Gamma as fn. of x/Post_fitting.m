% Post fitting_factor vector Simulation
DG = open('output/trial/dev1_sp/mat_files/file_0.3_DG.mat');

% Excluding the end points
x = DG.x(2:100)*0.0529;
v = DG.VDG(2:100)*19.2;

% Plot
figure(1);
plot(x,DG.fitting_factor);
ylabel('fitting factor');
xlabel('x(nm)');

figure(2);
plot(v,DG.fitting_factor);
ylabel('fitting factor');
xlabel('Potential(Volt)');
hold all;

% Polynomial Curve Fitting
coeff = polyfit(v,DG.fitting_factor,5);
disp(coeff);

yvalues = coeff(1).*(v.^5) + coeff(2).*(v.^4) + coeff(3).*(v.^3) + coeff(4).*(v.^2) + coeff(5).*v + coeff(6) ;

% Plot Comparison
figure(2);
plot(v,yvalues);
legend('Actual','fitted');