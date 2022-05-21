%%  Load raw data
clear all;
load data;

%%  Visualize the raw data 
L = length(MEP_raw(1,:));
t = (1:L)*ts - 0.2;     % make stimulus at t = 0;

figure(1);
plot(t*1000,MEP_raw');
title('All MEPs, full length');
set(gcf,'position',[200,200,500,200]);
xlabel('Time (ms)');
ylabel('MEP (\mu{V})');

%%  Obtain MUAP for this subject
post_range = 1100:1220;     % active MEP range
pre_range   = 1:950;        % background, for baseline noise evaluation
MEP = MEP_raw(:,post_range);
background = MEP_raw(:,pre_range);

v_init = mean(MEP(:,11:90));    % initilization of MUAP
%   the initialized MUAP is even shorter than the post MEP range, allowing
%   room for sliding in t-axis.
Nitr = 7;

[MUAP,record_MSE] = get_MUAP(MEP, v_init, Nitr);
figure(2);
	plot(record_MSE,'sk-');
	title('MUAP-learning process');
	xlabel('Iteration');
	ylabel('MSE of reconsruction');
	set(gcf,'position',[200,200,500,200]);

figure(3);hold on;
	plot((1:80)*0.2,MUAP,'k');
	title('Learned MUAP');
	xlabel('Time (ms)');
	set(gcf,'position',[200,200,500,200]);

%%  Decomposition & reconstruction
[~,MEP_est,Vpp_est] = get_decomp(MEP, MUAP);
Vpp_raw = range(MEP,2);

[~,idx_sort] = sort(STI);

%   plot IO curve
unique_STI = unique(STI);
M_Vpp_est = zeros(1,length(unique_STI));
M_Vpp_raw = zeros(1,length(unique_STI));
for i = 1 : length(unique_STI)
    sel_STI = (STI == unique_STI(i)); 
    M_Vpp_est(i) = median(Vpp_est(sel_STI));
    M_Vpp_raw(i) = median(Vpp_raw(sel_STI));
end

figure(4);
subplot(1,2,1); hold on;
	hold on;box on;grid on;
	scatter(STI(idx_sort),Vpp_est(idx_sort),'ko','filled');
	scatter(unique_STI,M_Vpp_est,'rs','filled');
	title('Vpp from resconstructed')
	xlabel('Stimulation strength (%)');
	set(gca,'yscale','log');
	ylabel('Vpp (est) (\mu{V})')
	ylim([1e-1,1e4]);

subplot(1,2,2)
	hold on;box on;grid on;
	scatter(STI(idx_sort),Vpp_raw(idx_sort),'ko','filled');
	scatter(unique_STI,M_Vpp_raw,'rs','filled');
	title('Vpp from raw MEP')
	xlabel('Stimulation strength (%)');
	set(gca,'yscale','log');
	ylabel('Vpp (est) (\mu{V})')
	ylim([1e-1,1e4]);
	set(gcf,'position',[100,100,800,350]);

