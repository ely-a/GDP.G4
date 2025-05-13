function [cdstr, utstr] = jd2str(jdtdb)

% convert Julian date to string equivalent
% calendar date and universal time

% input

%  jdtdb = Julian date

% output

%  cdstr = calendar date string
%  utstr = universal time string

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[month, day, year] = gdate(jdtdb);

% serial date number

sdn = datenum(year, month, day);

% create calendar date string

cdstr = datestr(sdn, 1);

% create universal time string

utstr = datestr(day - fix(day), 'HH:MM:SS.FFF');
