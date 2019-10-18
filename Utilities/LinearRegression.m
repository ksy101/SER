function [dl_Gd,dl_I,Gd_mean,Gd_std,I_mean,I_std] = LinearRegression(img_gd, img_I, Gd_id,I_id)

% img_gd = XL2H_hat_D(1:320,:,4); % image of Gd
% img_I = XL2H_hat_D(321:end,:,3); % image of Gd
idx4gd = [7 8 9  10 11 1 2 3 4]; % form 0.1 to 14 mg/cc
idx4i  = [8 9 10 11 12 1 3 4 5]; % form 0.1 to 14 mg/cc
idxfull = [5 6 7 8 9 10 11 1 2 3 4];
xx= [0.0001 0.0002 0.0005 0.001 0.002 0.005 0.010 0.012 0.015]*1e3; % x- ax
Gd_mean = zeros(1, length(idx4gd)); % mean
Gd_std  = zeros(1, length(idx4gd)); % std
I_mean = zeros(1, length(idx4gd)); % mean
I_std  = zeros(1, length(idx4gd)); % std
for i = 1:length(idx4gd)
    Gd_mean(i)  = mean(img_gd(Gd_id==idx4gd(i)))*1e3;
    Gd_std(i)   = std(img_gd(Gd_id==idx4gd(i)))*1e3;
    I_mean(i)  = mean(img_I(I_id==idx4i(i)))*1e3;
    I_std(i)   = std(img_I(I_id==idx4i(i)))*1e3; 
end
Gd_nonz = zeros(1,length(idxfull)); % pixel number in each roi
for i = 1:length(idxfull)
    img_mid     = img_gd(Gd_id==idxfull(i));
    img_mid(img_mid~=0)=1;   
    Gd_nonz(i)  = sum(img_mid); 
end   

% xx = xx(1:5);
% Gd_mean = Gd_mean(1:5);
% I_mean = I_mean(1:5);

dl_Gd = fitlm(xx,Gd_mean);
dl_I = fitlm(xx,I_mean);  

figure;
plot(xx,xx,'LineWidth',2)
grid on
hold on;plot(xx,Gd_mean,'o','LineWidth',2)
plot(xx,dl_Gd.Fitted,':','LineWidth',2)
plot(xx,I_mean,'*','LineWidth',2)
plot(xx,dl_I.Fitted,':','LineWidth',2)
ylabel('Measured concentration(mg/ml)')
xlabel('True concentration (mg/ml)')
set(gca,'xtick',[-1:1:15],'ytick',[-1:1:15])
axis([-1 15 -1 15]);
hl = legend('theoretical concentration','measured concentration - Gd','fitted measured concentration - Gd',...
    'measured concentration - I','fitted measured concentration - I');
set(hl,'Box','off');
% title('Decomposed Gd and I')