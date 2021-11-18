load('problem4.mat');



%%








%%



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

% 15:35 stopped working.. nause
% 16:00 continue working
% 17:45 stop working still in the train

% mu_max * m = 1/m * sum r^2

for i = 1:numel(m_list)

    for n = 1:numel(c_list)

        m = m_list(i);
        c = c_list(n);

        step = (2*c*mu_max)/m;

 

            [y, e, w] = doLMSFiltering(m,step,chestECG,abdomenECG);

            MSE = evaluateResult(e);

            if MSE < best_mse
                best_mse = MSE;
                best_m = m;
                best_c = c;
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
         lms = dsp.LMSFilter(m, 'StepSizeSource', 'Input Port');
         [y,e,w] = lms(r,x,step);
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

    mse = immse(fetalECG(INITIAL_REJECTION:end), v(INITIAL_REJECTION:end));
    end
end
