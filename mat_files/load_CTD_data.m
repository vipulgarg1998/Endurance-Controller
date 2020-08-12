function [depth, density_data] = load_CTD_data(filename, path)
    %%%%%%%%%%%%%%%%% read data from CTD dataset %%%%%%%%%%%%%%%%%%%%%
    data = importdata(filename) ;
    CTD_data = data.data(:,:) ;
    CTD_data_size = size(CTD_data,1) ;
    depth = CTD_data(:,1) ;

    density_data = zeros(CTD_data_size,1) ;
    current_dir = cd(path) ;
    files = dir('*.m') ;
    names = cell(size(files)) ;
    parfor i = 1:max(size(files))
        names{i} = files(i).name ;
    end
    poolobj = gcp;
    addAttachedFiles(poolobj,names') ;
    cd(current_dir) ;
    parfor i = 1:max(size(CTD_data))
        % density = f( salinity, temperature, pressure)
        density_data(i) = sw_dens(CTD_data(i,5), CTD_data(i,3), CTD_data(i,2)) ;
    end
end