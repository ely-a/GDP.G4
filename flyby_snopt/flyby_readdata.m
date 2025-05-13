function [fid, jdtdbi1, ddays1, jdtdbi2, ddays2, jdtdbi3, ddays3, ...
    ip1, ip2, ip3, otype] = flyby_readdata(filename)

% read flyby simulation definition data file

% NOTE: all angular elements are returned in radians

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fbalt_lwr fbalt_upr rkcoef

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)

    clc; home;

    fprintf('\n\n error: cannot find this file!!');

    pause

    return;

end

% read 96 lines of data file

for i = 1:1:96

    cline = fgetl(fid);

    switch i

        case 18

            % type of trajectory optimization

            otype = str2double(cline);

        case 24

            % departure calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jdtdbi1 = julian(month, day, year);

        case 27

            % search boundary for departure date (days)

            ddays1 = str2double(cline);

        case 42

            % departure planet

            ip1 = str2double(cline);

        case 48

            % flyby calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jdtdbi2 = julian(month, day, year);

        case 51

            % search boundary for flyby date (days)

            ddays2 = str2double(cline);

        case 66

            % flyby planet

            ip2 = str2double(cline);

        case 69

            % lower bound for flyby altitude (kilometers)

            fbalt_lwr = str2double(cline);

        case 72

            % upper bound for flyby altitude (kilometers)

            fbalt_upr = str2double(cline);

        case 78

            % arrival calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jdtdbi3 = julian(month, day, year);

        case 81

            % search boundary for arrival date (days)

            ddays3 = str2double(cline);

        case 96

            % arrival celestial body

            ip3 = str2double(cline);

    end

end

% initialize rkf78 method

rkcoef = 1;

fclose(fid);

