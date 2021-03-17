function [X] = create_labels(data)
% [X] = createLabels(data)
%   Creates labels using the information provided in the data.
num_elements = size(data.station_counts,1);
num_sensors  = size(data.station_counts,2);
X = zeros(num_elements, num_sensors);
for i=1:length(data.station_ids)
    ind_start = round(minutes(data.inc_timestamp(i,:)-...
        datetime(data.station_times(1,:)))/5);
    ind_end = round(data.inc_duration(i)/5)+ind_start;
    [~,ind_station,~] = intersect(data.station_ids,data.inc_station_id(i));
    X(ind_start:ind_end, ind_station) = 1;
end
end