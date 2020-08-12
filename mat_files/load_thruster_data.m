function [thruster_data, power_data] = load_thruster_data(path)
    data = xlsread(path) ;
    % lbf to newton conversion factor
    lbf_to_N = 4.44822 ;
    % extract data from file
    thruster_data = data(:,1).*lbf_to_N ;   % in Newtons
    power_data = data(:,2) ;              % in Watts

    [~, ind] = unique(thruster_data, 'rows');
    % duplicate indices
    duplicate_ind = setdiff(1:size(thruster_data, 1), ind);
    % remove duplicate values
    thruster_data(duplicate_ind) = [] ;
    power_data(duplicate_ind) = [] ;

end