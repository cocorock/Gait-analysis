
% Functions_rev/apply_butterworth_filter.m
function filtered_data = apply_butterworth_filter(data, sampling_freq, cutoff_freq)
    % APPLY_BUTTERWORTH_FILTER Applies a zero-lag Butterworth low-pass filter.
    %   filtered_data = APPLY_BUTTERWORTH_FILTER(data, sampling_freq, cutoff_freq)
    %   applies a 4th order Butterworth low-pass filter to the input data.
    %   filtfilt is used to ensure zero phase distortion.
    %
    %   Inputs:
    %     data: Input data (can be a vector or a matrix where filtering is
    %           applied column-wise).
    %     sampling_freq: Sampling frequency of the data (Hz).
    %     cutoff_freq: Cutoff frequency for the low-pass filter (Hz).
    %
    %   Output:
    %     filtered_data: The filtered data.

    if nargin < 3
        error('Not enough input arguments. Usage: apply_butterworth_filter(data, sampling_freq, cutoff_freq)');
    end

    if cutoff_freq >= sampling_freq / 2
        warning('Cutoff frequency is at or above Nyquist frequency. No filtering applied.');
        filtered_data = data;
        return;
    end

    % Design a 4th order Butterworth filter
    order = 4;
    [b, a] = butter(order, cutoff_freq / (sampling_freq / 2), 'low');

    % Apply the filter using filtfilt for zero phase distortion
    filtered_data = filtfilt(b, a, data);
end
