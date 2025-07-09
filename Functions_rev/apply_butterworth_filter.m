% Functions_rev/apply_butterworth_filter.m
function filtered_data = apply_butterworth_filter(data, sampling_freq, cutoff_freq)
    % APPLY_BUTTERWORTH_FILTER Applies a zero-lag Butterworth low-pass filter.
    %   filtered_data = APPLY_BUTTERWORTH_FILTER(data, sampling_freq, cutoff_freq)
    %   applies a 4th order Butterworth low-pass filter to the input data.
    %   filtfilt is used to ensure zero phase distortion.
    %   This version includes circular padding to reduce edge effects.
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

    % --- Padding to reduce edge effects ---
    pad_len = 60; % Length of padding on each side

    % Ensure data has enough length for padding
    if size(data, 1) < pad_len * 2 + 1
        warning('Data length is too short for specified padding. Skipping padding.');
        padded_data = data;
    else
        % Circular padding for multi-column data
        % Take last 'pad_len' elements, original data, first 'pad_len' elements
        padded_data = [data(end-pad_len+1:end, :); data; data(1:pad_len, :)];
    end

    % Apply the filter using filtfilt for zero phase distortion
    filtered_padded = filtfilt(b, a, padded_data);

    % Extract the original portion of the filtered data
    if size(data, 1) < pad_len * 2 + 1
        filtered_data = filtered_padded; % No padding was applied, return as is
    else
        filtered_data = filtered_padded(pad_len+1 : pad_len+size(data, 1), :);
    end
end