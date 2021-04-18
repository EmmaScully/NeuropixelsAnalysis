function average_matrix = average_flashes(channels, maxflashrate, data, threshold)

average_matrix = zeros(channels, maxflashrate);
number_of_flashes = ceil(length(data)/maxflashrate);
blank_matrix = zeros(maxflashrate, number_of_flashes);
b = 1;
for c = 1:channels
    for d = 1:length(data)
        for a=1:maxflashrate
            blank_matrix(b,a) = data(channels, d);
            if data(channels, d) < threshold
                b = b+1;                
            end
        end
    end
    average_matrix(c,:) = mean(blank_matrix);
end





end