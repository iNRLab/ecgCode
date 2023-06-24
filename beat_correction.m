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

% Below, we first detect ectopic beats based on the presence of NPN or PNP patterns in the dRR series 
% and replace the central value with the average of its surrounding values. 
% Then we detect missed and extra beats by comparing the current RR interval or the 
% sum of two successive RR intervals with the median of the surrounding RR intervals. 
% For a missed beat, the RR interval is replaced with the median value. For an extra beat, 
% two successive RR intervals are each replaced with half of the median value.


function corrected_data = automatic_correction(raw_data)
    window_size = 90;
    med_window_size = 10;  % window for the median calculation
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

        % detect ectopic beats
        if i >= 2 && i <= length(dRR) - 1
            if ((dRR(i-1) < 0 && dRR(i) > 0 && dRR(i+1) < 0) || ...  % NPN pattern
                (dRR(i-1) > 0 && dRR(i) < 0 && dRR(i+1) > 0)) && ...  % PNP pattern
                abs(dRR(i)) > Th
                corrected_data(i) = (raw_data(i-1) + raw_data(i+1)) / 2;
            end
        end

        % detect missed or extra beats
        if i > med_window_size/2 && i <= length(dRR) - med_window_size/2
            medRR = median(raw_data((i - med_window_size/2):(i + med_window_size/2)));

            if raw_data(i) > 1.5 * medRR  % condition for a missed beat
                corrected_data(i) = medRR;
            elseif i < length(raw_data) && ...  % ensure that we are not at the last data point
                   raw_data(i) + raw_data(i+1) < 0.8 * medRR  % condition for an extra beat
                corrected_data(i) = medRR / 2;
                corrected_data(i+1) = medRR / 2;
            end
        end
    end
end

