function target = DBMerge(source,target)

%DBMerge - Merge databases.
%
%  USAGE
%
%    target = DBMerge(source,target)
%
%    source         database where the data is read
%    target         database where the data is added
%
%  SEE
%
%    See also DBCreate, DBDuplicate.
%

% Copyright (C) 2007-2016 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end
if ~isastring(source),
  error('Incorrect source database name (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end
if ~isastring(target),
  error('Incorrect target database name (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end

% Copy database contents
try
	h = mym(['insert into ' target '.' 'figures select * from ' source '.figures']);
	h = mym(['insert into ' target '.' 'variables select * from ' source '.variables']);
catch
   error('FMAToolbox:DBMerge:mergeDB',['Could not merge databases ''' source ''' and ''' target '''.']);
end
