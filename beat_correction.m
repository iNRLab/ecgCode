% Threshold based beat correction algorithm
function corrected_data = threshold_based_correction(raw_data, threshold)
    window_size = 90;  % assuming a window of 90 RR intervals
    corrected_data = raw_data;
    local_average = medfilt1(raw_data, window_size);

    for i = 1:length(raw_data)
        if abs(raw_data(i) - local_average(i)) > threshold
            % replace the value with the average of the surrounding values
            % in the case where the artifact is at the beginning or end, replace with the closest valid value
            if i == 1
                corrected_data(i) = raw_data(i+1);
            elseif i == length(raw_data)
                corrected_data(i) = raw_data(i-1);
            else
                corrected_data(i) = (raw_data(i-1) + raw_data(i+1)) / 2;
            end
        end
    end
end

% Automatic beat correction algorithm
function corrected_data = automatic_correction(raw_data)
    window_size = 90;
    corrected_data = raw_data;
    dRR = diff(raw_data);  % dRR series

    for i = 1:length(dRR)
        % calculate quartile deviation of surrounding 90 beats and multiply by 5.2
        if i < window_size/2
            Th = 5.2 * iqr(dRR(1:(i + window_size/2)));
        elseif i > length(dRR) - window_size/2
            Th = 5.2 * iqr(dRR((i - window_size/2):end));
        else
            Th = 5.2 * iqr(dRR((i - window_size/2):(i + window_size/2)));
        end

        % identify the beats that exceed the threshold
        if abs(dRR(i)) > Th
            if i == 1
                corrected_data(i) = raw_data(i+1);
            elseif i == length(raw_data)
                corrected_data(i) = raw_data(i-1);
            else
                corrected_data(i) = (raw_data(i-1) + raw_data(i+1)) / 2;
            end
        end
    end
end
