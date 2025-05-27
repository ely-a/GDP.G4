function load_moon_data()
    folder = 'Ephemerides';
    files = dir(fullfile(folder, '*.txt'));
    orbits = {};
    velocities = {};
    moon_names = {};
    times_all = {};
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
        pos_data = [];
        vel_data = [];
        times = [];
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, '$$EOE')
                break;
            end
            parts = split(line, ',');
            if numel(parts) >= 8
                xyz = str2double(parts(3:5));
                vxyz = str2double(parts(6:8));
                t = str2double(strtrim(parts{1})); % Assume timestamp is in JD
                if all(~isnan(xyz)) && all(~isnan(vxyz)) && ~isnan(t)
                    pos_data = [pos_data; xyz']; %#ok<AGROW>
                    vel_data = [vel_data; vxyz']; %#ok<AGROW>
                    times = [times; t]; %#ok<AGROW>
                end
            end
        end
        fclose(fid);

        if ~isempty(pos_data)
            orbits{end+1} = pos_data;
            velocities{end+1} = vel_data;
            moon_names{end+1} = erase(files(k).name, '.txt');
            times_all{end+1} = times;
            max_len = max(max_len, size(pos_data,1));
        end
    end

    % Pad and stack
    for i = 1:length(orbits)
        n = size(orbits{i}, 1);
        if n < max_len
            orbits{i}(end+1:max_len, :) = repmat(orbits{i}(end, :), max_len - n, 1);
            velocities{i}(end+1:max_len, :) = repmat(velocities{i}(end, :), max_len - n, 1);
            times_all{i}(end+1:max_len, 1) = repmat(times_all{i}(end), max_len - n, 1);
        end
    end

    orbits_mat = cat(3, orbits{:});
    velocities_mat = cat(3, velocities{:});
    save('moon_orbits.mat', 'orbits_mat', 'velocities_mat', 'moon_names', 'max_len', 'times_all');
    disp('Saved moon_orbits.mat');
end
