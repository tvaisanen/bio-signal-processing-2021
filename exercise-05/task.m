% The sampling rate is 2000 Hz 
FS = 2000;

% Load the signals from data.mat into the struct 'data'
% << insert loading code here >>
data = load('data.mat').data;
% Number of segments
N = numel(data);  
idx = 1:5;

% Calculate average force of each segment (1xN vector)
AF = arrayfun(@(i) mean(data(i).force), idx);

% Calculate EMG dynamic range in each segment (1xN vector)
DR = arrayfun(@(i) max(data(i).EMG) - min(data(i).EMG), idx);

% % Calculate EMG mean squared value in each segment (1xN vector)
MS = arrayfun(@(i) sum(data(i).EMG .^ 2)/numel(data(i).EMG), idx);
% 
% 
r = 1:200;
y = data(2).EMG;
%s = sign(y)
% c = s~=[s(2:end);s(end)];


% crossings = find(c(r));
% figure(1)
% clf
% plot(data(1).t(r),y(r));
% %hold on
% %plot(r,0.1+s*0.2)
% hold on
% plot(data(1).t(crossings),zeros(1,numel(crossings)), 'o')

zcr_2 = @(s) s~=[s(2:end);s(end)];
zcr_1 = @(y) zcr_2(sign(y));
zcr_0 = @(i) numel(find(zcr_1(data(i).EMG)))/(numel(data(i).t)/FS);
ZCR = arrayfun(zcr_0, idx);


% get_diffs = @(c) sum(c(c~=[c(2:end);2]))/numel(c);
% get_signs = @(i) get_diffs(sign(diff(data(i).EMG)));
% 
% % % Calculate EMG zero crossing rate in each segment (1xN vector)
% TCR = arrayfun(get_signs, idx);

%% TCR

TCR = zeros(1,N);


for i = 1:N
    
    y = data(i).EMG;
    s = sign(diff(y));
    s1 = [s;s(end)];
    s2 = [s(1);s];

    ve_diff = (s2 - s1)  ~= 0;   
    ve_mult = (s1 .* s2) ~= 1;
    ve_turn = ve_diff & ve_mult;
    ve_coord= find(ve_turn);
    
    tcr_coord = find(abs(diff(y(ve_coord))) > 0.1);
    
    TCR(i) = FS*numel(tcr_coord)/(data(i).length);
    
    if i == 1
        clf
        plot(data(i).t,y)
        hold on
        plot(data(i).t(tcr_coord),y(tcr_coord),'r*');
    end
    
end

% % Calculate the linear model coefficients for each parameter
% % The models are in the form: parameter(force) = constant + slope * force,
% % and the coefficients are stored in a 1x2 vectors: p_<param> = [slope constant]
% % For example, p_DR(1) is the slope and p_DR(2) is the constant of the linear model mapping the average force into the dynamic range
% % You can use the 'polyfit' (or the 'regress') command(s) to find the model coefficients
p_DR  = polyfit(DR,idx,1);
p_MS  = polyfit(MS,idx,1);
p_ZCR = polyfit(ZCR,idx,1);
p_TCR = polyfit(TCR,idx,1);

%% CORRELATIONS
% % Calculate correlation coefficients between the average forces and each of the parameters using 'corr'
c_DR = corr(AF',DR');
c_MS = corr(AF',MS');
c_ZCR = corr(AF',ZCR');
c_TCR = corr(AF',TCR');