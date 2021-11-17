% Required reading
% Read chapters 3.6-3.9 -- especially the chapter 3.9 -- from the course book.

% Assignment Description

% Your task is to implement an adaptive filtering solution to obtain fetal 
% ECG by cancelling out the maternal ECG from the measurement.

% Data
% The file 'signals.mat' contains all the signals of the assignment.
% The variables contained in the file are listed below:
% Fetus signal: 'fhb' (pure fetus signal that is used to assess analysis results)
% Mother’s chest signal: 'mhb' (pure mother’s signal that is used to assess analysis results)
% Abdomen signals 'abd_sig1', 'abd_sig2' and 'abd_sig3' (mixed from fetus and mother signals)
% A real respiration movement signal and R-to-R interval sequence from ECG signal: 
% 'RespReference' and 'RRiInput'
% The sampling frequency of 'fhb', 'mhb', 'abd_sig1', 'abd_sig2' and 'abd_sig3'
% is 1000 Hz, and they are stored as 20000x1 vectors. Additionally, the sampling 
% rate of 'RespReference' and 'RRiInput' is 4 Hz, and they are 1596x1 vectors.

% The sampling rates are 1000 Hz
FS = 1000;

%% Task 1: Pure Sum of Signals

% Load problem1.mat to have access to variables abd_sig1 and mhb

load("problem1.mat")

% Calculate sample timing vector in seconds starting from 0
t = 0:1/FS:20-1/FS;

% Estimate c2 from abd_sig1 and mhb using the scalar projection formula
c2_projection = dot(abd_sig1, mhb)/dot(mhb,mhb);


% Estimate c2 from abd_sig1 and mhb using the pseudoinverse function (pinv)
c2_pinv = pinv(mhb) * abd_sig1;
% Estimate c2 from abd_sig1 and mhb using the backslash operator (\)
c2_operator = mhb\abd_sig1;

% Calculate the fetal ECG by cancelling out the scaled mother's ECG using projection based estimation coefficient
fetus = abd_sig1 - c2_projection * mhb;


figure

subplot(311)
plot( t, abd_sig1, 'b' )
hold on
plot( t, mhb, 'r--' )
legend('abd\_sig1 (x = x_1 + c_2 x_2)','mhb (y = x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(312)
plot( t, c2_projection * mhb, 'b'  )
legend('scaled mhb (c_2 x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(313)
plot( t, fetus, 'b' )
hold on
plot( t, fhb, 'r--' );
legend('estimated fetus (x - y = x - c_2 x_2)', 'true fetus')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);


%% Task 2:


load('problem2.mat')
% Calculate sample timing vector in seconds starting from 0
t = 0:1/FS:20-1/FS;


% Estimate the time lag using cross correlation
% (Calculate cross correlation using the xcorr function and then
% use the max function to find the lag giving maximal correlation)
[m,d] = max(xcorr(mhb_ahead,abd_sig1));
d = numel(abd_sig1) - d;

% Shift the chest ECG mhb_ahead back in time d samples padding with nearest value

% not correct :: mhb = [mhb_ahead(end)*ones(1,d-1), mhb_ahead(d:end)']';
% not correct :: mhb = [mhb_ahead(1)  *ones(1,d),   mhb_ahead(d+1:end)']';
% not correct :: mhb = [mhb_ahead(1)  *ones(1,d-1), mhb_ahead(d:end)']';

% not correct :: mhb = [mhb_ahead(d:end)',          mhb_ahead(end)*ones(1,d-1)]';
% not correct :: mhb = [mhb_ahead(d+1:end)',        mhb_ahead(end)*ones(1,d)]';

% mhb = [mhb_ahead(d+1:end)', mhb_ahead(end)*ones(1,d)]';
% mhb = [mhb_ahead(end)*ones(d+1,1);mhb_ahead(d+2:end)];
% This got me right c2 value
mhb = [mhb_ahead(1)*ones(d,1);mhb_ahead(1:end-d)];

% Estimate c2 from abd_sig1 and mhb
c2 = dot(abd_sig1,mhb)/dot(mhb,mhb);

% Calculate the fetal ECG by cancelling out the scaled mother's ECG using projection based estimation coefficient
fetus = abd_sig1 - c2 * mhb;

N = 500
clf
plot(t(end-N:end),abd_sig1(end-N:end))
hold on
plot(t(end-N:end),mhb(end-N:end))
hold on

figure

subplot(311)
plot( t, abd_sig1, 'b' )
hold on
plot( t, mhb, 'r--' )
legend('abd\_sig1 (x = x_1 + c_2 x_2)','mhb (y = x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(312)
plot( t, c2 * mhb )
legend('scaled mhb (c_2 x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(313)
plot( t, fetus, 'b' )
hold on
plot( t, fhb, 'r--' )
legend('fetus (x - y = x - c_2 x_2)', 'true fetus')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

%% Task 3:

load('problem3.mat')

% Calculate sample timing vector in seconds starting from 0
t = (0:1/FS:20-1/FS);
t2 = (0:1/(FS*2):20-1/(2*FS));
% Estimate the time lag using cross correlation with the 'xcorr' function
% Fit a spline to the cross correlation using 'spline' function, and then 
% find the delay with maximum correlation using 'fnmin'
% NOTE: to use minimization function for maximization, please invert the objective function!
[r, lags] = xcorr(mhb_ahead,abd_sig1);
pp = spline(lags,-r);
[a,d] = fnmin(pp,[-5,1]);
d = -d;
% fnplt(pp,[-5,1])
% Shift the chest ECG mhb_ahead back in time d samples
% Use linear interpolation with extrapolation with the function 'interp1'
integer_offset = floor(d);
partial_offset = d - integer_offset; % get rid of the integer

mhb = [mhb_ahead(1)*ones(integer_offset,1);
       mhb_ahead(1:end-integer_offset)];

t2 = t + partial_offset*0.0001;
mhb = interp1(t, mhb, t2, 'linear', 'extrap')';
clf
N = 200;

subplot(211)
plot(t(1:N),abd_sig1(1:N))
hold on
%plot(t(1:N),circshift(mhb_ahead(1:N),3))

hold on
plot(t(1:N),mhb(1:N))
legend('interpolated','mhb')

subplot(212)

plot(t(end-N:end),abd_sig1(end-N:end))
hold on
%plot(t(1:N),circshift(mhb_ahead(1:N),3))

hold on
plot(t(end-N:end),mhb(end-N:end))
legend('interpolated','mhb')


% Estimate c2 from abd_sig1 and mhb
c2 = dot(abd_sig1,mhb)/dot(mhb,mhb);

% Calculate the fetal ECG by cancelling out the scaled mother's ECG using projection based estimation coefficient
fetus = abd_sig1 - c2 * mhb;

figure

subplot(311)
plot( t, abd_sig1, 'b' )
hold on
plot( t, mhb, 'r--' )
legend('abd\_sig1 (x = x_1 + c_2 x_2)','mhb (y = x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(312)
plot( t, c2 * mhb )
legend('scaled mhb (c_2 x_2)')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

subplot(313)
plot( t, fetus, 'b' )
hold on
plot( t, fhb, 'r--' )
legend('fetus (x - y = x - c_2 x_2)', 'true fetus')
xlabel('t [s]')
ylabel('amplitude [a.u.]')
ylim([-2 4]);

%% Task 4:

load('problem4.mat');

MU_MAX = 0.05;
FILTER_LENGHTS = [1 5 11 15 21 31 51 101]';
ADAPTATION_RATES = (0.1:0.1:1)';

[best_m, best_c, best_w, best_mse] = findBestFilterParameters(mhb_ahead_PI, abd_sig1, fhb, FILTER_LENGHTS, ADAPTATION_RATES, MU_MAX)

function [best_m, best_c, best_w, best_mse] = findBestFilterParameters(chestECG, abdomenECG, fetalECG, m_list, c_list, mu_max)
% This function finds the best LMS filter parameter combination from the
% given lists using two inner functions.
% To be completed by you!
%
% INPUTS:
%   chestECG    ECG from the chest, maternal ECG only, reference input, Nx1
%   abdomenECG  ECG from the abdomen, fetal and maternal mixed, primary input, Nx1
%   fetalECG    ECG from the fetus alone, signal of interest, Nx1 (cannot be measured directly, but is given for evaluation here)
%   m_list      list of filter lengths/orders to test, Mx1
%   c_list      list of step size fractions to test, Cx1 (each >0 & <1)
%   mu_max      maximum step size, scalar
%
% OUTPUTS:
%   best_m      the best filter length (from the m_list), scalar
%   best_c      the best fraction of mu_max (from the m_list), scalar
%   best_w      the best filter coefficients, best_m x 1 vector
%   best_mse    the lowest mean squared error obtained with the best parameters, scalar

% When evaluating the results in evaluateResult(), skip this many samples from the beginning to avoid initial adaptation transient
INITIAL_REJECTION = 2000;
best_m = 0;
best_c = 0;
best_w = 0;
best_mse = 10000;

% Here you go through all the possible combinations of filter lengths in m_list and learning rate fractions in c_list selecting the best performing one
% << INSERT YOUR CODE HERE >>
% 



for m = 1:numel(m_list)
    for n = 1:numel(c_list)
        L = m_list(m);
        C = c_list(n);

        step = 2*C/sum(power(c_list,2));
        
        [y, e,w] = doLMSFiltering(L,step,chestECG,abdomenECG);
        error = evaluateResult(y)

        if error < best_mse
            best_mse = error;
            best_m = L;
            best_c = C;
            best_w = w;
        end

    end
end

    function [y,e,w] = doLMSFiltering(m,step,r,x)
    % Does the actual LMS filtering.
    % To be completed by you!
    %
    % INPUTS:
    %   m       filter length
    %   step    LMS learning rule step size
    %   r       reference input (to be filtered)
    %   x       primary observed signal
    %
    % OUTPUTS:
    %   y       filtered signal r
    %   e       filter output, estimate of the signal of interest v

    % Create the dsp.LMSFilter object and use it to filter the input data
    % << INSERT YOUR CODE HERE >>

         lms = dsp.LMSFilter(m, 'StepSizeSource', 'Input Port');
         [y,e,w] = lms(x,r,step);
     end

    function mse = evaluateResult(v)
    % Calculates the mean squared error between the filtered signal v and
    % the known fetal ECG.
    %
    % NOTE1:    Skip INITIAL_REJECTION number of samples in the beginning of both signals to not include initial adaptation transient
    % 
    % NOTE2:    This nested function can access the desired output value in fetalECG directly!
    %
    % INPUTS:
    %   v       estimate of the signal of interest 
    %
    % OUTPUTS:
    %   mse     mean squared error between v and fetalECG    

    % You can call the 'immse' function for the signals without the initial rejection parts
    % << INSERT YOUR CODE HERE >>
    mse = immse(fetalECG(INITIAL_REJECTION+1:end), y(INITIAL_REJECTION+1:end))
    end
end








