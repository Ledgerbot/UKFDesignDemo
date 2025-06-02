function varargout = UKFDesignDemoVer
% UKFDESIGNDEMOVER displays the UKF Design Demo information.
%   UKFDESIGNDEMOVER displays the information to the command prompt.
%
%   A = UKFDESIGNDEMOVER returns in A the sorted struct array of  
%   version information for the UKF Design Demo.
%     The definition of struct A is:
%             A.Name      : toolbox name
%             A.Version   : toolbox version number
%             A.Release   : toolbox release string
%             A.Date      : toolbox release date
%
%   M. Kutzer 02Jun2025, USNA

% Updates


A.Name = 'UKF Design Demo';
A.Version = '1.0.0';
A.Release = '(R2020b)';
A.Date = '02-Jun-2025';
A.URLVer = 1;

msg{1} = sprintf('MATLAB %s Version: %s %s',A.Name, A.Version, A.Release);
msg{2} = sprintf('Release Date: %s',A.Date);

n = 0;
for i = 1:numel(msg)
    n = max( [n,numel(msg{i})] );
end

fprintf('%s\n',repmat('-',1,n));
for i = 1:numel(msg)
    fprintf('%s\n',msg{i});
end
fprintf('%s\n',repmat('-',1,n));

if nargout == 1
    varargout{1} = A;
end