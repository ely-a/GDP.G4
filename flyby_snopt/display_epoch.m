function display_epoch(jdate, jtype)

% display calendar date and time

% input

%  jdate = julian day

%  jtype = type of julian day (1 = TDB, 2 = UTC)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cdate, ctime] = jd2str(jdate);

switch jtype

    case 1

        fprintf('\nTDB calendar date     ');

        disp(cdate);

        fprintf('TDB time              ');

        disp(ctime);

        fprintf('TDB Julian day       %16.7f\n', jdate);

    case 2

        fprintf('\nUTC calendar date     ');

        disp(cdate);

        fprintf('UTC time              ');

        disp(ctime);

        fprintf('UTC Julian day       %16.7f\n', jdate);

end
