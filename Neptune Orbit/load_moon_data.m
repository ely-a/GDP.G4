function load_moon_data()
    folder = 'Ephemerides';
    files = dir(fullfile(folder, '*.txt'));
    orbits = {};
    moon_names = {};
    max_len = 0;

    for k = 1:length(files)
        fname = fullfile(folder, files(k).name);
        fid = fopen(fname, 'r');

        % Skip to $$SOE
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, '$$SOE')
                break;
            end
        end

        % Read ephemeris lines
        data = [];
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, '$$EOE')
                break;
            end
            parts = split(line, ',');
            if numel(parts) >= 5
                xyz = str2double(parts(3:5));
                if all(~isnan(xyz))
                    data = [data; xyz']; %#ok<AGROW>
                end
            end
        end
        fclose(fid);

        if ~isempty(data)
            orbits{end+1} = data;
            moon_names{end+1} = erase(files(k).name, '.txt');
            max_len = max(max_len, size(data,1));
        end
    end

    % Pad and stack
    for i = 1:length(orbits)
        n = size(orbits{i}, 1);
        if n < max_len
            pad = repmat(orbits{i}(end, :), max_len - n, 1);
            orbits{i} = [orbits{i}; pad];
        end
    end

    orbits_mat = cat(3, orbits{:});
    save('moon_orbits.mat', 'orbits_mat', 'moon_names', 'max_len');
    disp('Saved moon_orbits.mat');
end