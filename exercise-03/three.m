% Load problem3.mat to have access to variables abd_sig1 and mhb_ahead

load('problem3.mat');
% The sampling rates are 1000 Hz
FS = 1000;

% Calculate sample timing vector in seconds starting from 0
t =  (0:1/FS:20-1/FS);

% Estimate the time lag using cross correlation with the 'xcorr' function
% Fit a spline to the cross correlation using 'spline' function, and then find the delay with maximum correlation using 'fnmin'
% NOTE: to use minimization function for maximization, please invert the objective function!
[r, lags] = xcorr(mhb_ahead,abd_sig1);
pp = spline(lags,-r);
[a,d] = fnmin(pp,[-100,100])
d = -d

% Shift the chest ECG mhb_ahead back in time d samples
% Use linear interpolation with extrapolation with the function 'interp1'
t2 = t - d*(1/FS);

mhb = interp1(t,mhb_ahead, t2, 'linear','extrap')';
% Estimate c2 from abd_sig1 and mhb
c2 = dot(abd_sig1,mhb)/dot(mhb,mhb)


fetus = abd_sig1 - c2 * mhb;