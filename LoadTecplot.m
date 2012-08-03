function [dat] = LoadTecplot(fname,varargin)
% Load Tecplot ASCII FEBLOCK file into a structure.
% 
%  [dat] = LoadTecplot(fname,{options})
%
% Inputs:
%     fname : Filename to load (if empty or not found then a dialog will be presented)
%     options : string options
%               'plot' - plot surface data
%               'save' - save the data in Matlab format
%               'delete' - delete the source file, will force a save in Matlab format
%               'quiet' - don't display any text info while loading
%
% Outputs:
%     dat : Tecplot dat in a block structure of variables
%           each zone is a field of this and each variable is a field of each zone.
%           Additional fields include;
%           source : filename of source data
%           title : title field from Tecplot file
%           vars : variables read in
%           zones : list of zones 
%           min, max, mean, std : replicated structure with info about each zone variable
%
% Currently only loads surface data; TRI and QUAD.
%
% If no output variable is given (nargout=0) then save and plot are automatically activated.
%
% Examples:
%
%    % Open a GUI to select a file and plot the contents
%    LoadTecplot
%
%    % Load a file, save the Matlab structure and plot the data
%    dat = LoadTecplot('data.dat','plot','save');
%



% file to load
if (~exist('fname','var') || ~exist(fname,'file')),
  [filename, pathname] = uigetfile('*.dat', 'Pick an Tecplot-file');
  fname = [pathname filename];
end;

% defaults
if (nargout==0),  % Save and plot automatically if no output variable is specified
  saveem = true;
  plotem = true;
else
  saveem = false;
  plotem = false;
end;
deleteem = false;
quiet = 'false';
for ii = 1:length(varargin),
  switch(lower(varargin{ii})),
  case('plot')
    plotem = true;
  case('save')
    saveem = true;
  case('delete')
    deleteem = true;
    saveem = true;
  case('del')
    deleteem = true;
    saveem = true;
  case('quiet')
    quiet = true;
  otherwise
    error([mfilename ':input'],'Unknown Option %s',varargin{ii});
  end;
end;

% open Tecplot file and load
fid = fopen(fname,'r');
if (fid<=0), error([mfilename ':io'],'Failed to load %s',fname); end;
dat.source = fname;

% Get title
%TITLE = "title"
ll = fgetl(fid);
res = regexp(ll,'TITLE *= *"(?<title>.*)"', 'names');
if ~isempty(res),
  dtitle = res(1).title;
  if (~quiet), disp(sprintf('Title:%s',dtitle)); end;
end;
dat.title = 'lala';

% Get variable names
%VARIABLES  = X, Y, Z, "pressure", "density", "y-plus", "x-wall-shear", "y-wall-shear", "z-wall-shear", "x-face-area", "y-face-area", "z-face-area", "cell-wall-distance"
ll = fgetl(fid);
res = regexp(ll,'VARIABLES *= *(?<vars>.*) *', 'names');
if ~isempty(res),
  res = regexp([res(1).vars ','],'[" ]*(?<name>[^,]*?)[" ]*,', 'tokens');
  for ii = 1:length(res),
    oname{ii} = res{ii}{1};
    vname{ii} = regexprep(res{ii}{1},'[\W]+','_');
    %dat.(vname{ii}) = struct;
  end;
  dat.vars = vname;
end;

% plot setup
if (plotem)
  % Look for pressure variable
  IV = strmatch('pressure',vname,'exact');
  if isempty(IV),
    IV = min(4,length(vname));  % Use first variable
  end;
  % Open figure
  figure;
  axis('equal');
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title(vname{IV});
  lighting('flat');
  hold on;
  cameratoolbar('show','setmode','orbit');
end;

% Get each Zone
%ZONE T="aft/a_exit", N=262, E=431, ET=TRIANGLE, F=FEBLOCK
ll = fgetl(fid);
count = 0;
while ischar(ll),
  % Read Zone header
  res = regexp(ll,'ZONE +T="(?<zone>[^"]+)", N=(?<N>[^"]+), E=(?<E>[^"]+), ET=(?<ET>[^"]+), F=(?<F>[^"]+)', 'names');
  if ~isempty(res),
    count = count + 1;

    % Check type
    if ~strcmpi(res(1).F,'FEBLOCK'), error([mfilename ':format'],'Unknown structure %s',res(1).F); end;

    % Zone info
    zone = regexprep(res(1).zone,'[\W]+','_');
    dat.zones{count} = zone;
    if (~quiet), disp(sprintf('\t%3i: %s %s nodes %s %s',count,res(1).zone,res(1).N,res(1).E,res(1).ET)); end;
    dat.(zone) = struct;
    dat.(zone).name = res(1).zone;
    dat.(zone).N = str2double(res(1).N);
    dat.(zone).E = str2double(res(1).E);
    dat.(zone).ET= res(1).ET;
    dat.(zone).F = res(1).F;

    % read data blocks
    for iv = 1:length(vname),
      tmp  = fscanf(fid,'%f',dat.(zone).N);
      dat.(zone).max.(vname{iv}) = max(tmp(:));
      dat.(zone).min.(vname{iv}) = min(tmp(:));
      dat.(zone).mean.(vname{iv})= mean(tmp(:));
      dat.(zone).std.(vname{iv}) = std(tmp(:));
      dat.(zone).(vname{iv}) = tmp; 
    end;

    switch(dat.(zone).ET),
    case('TRIANGLE'),
      dat.(zone).tri = reshape(fscanf(fid,'%f',dat.(zone).E*3),[3,dat.(zone).E])';
      % create list of edges
      edge = reshape(dat.(zone).tri(:,[1 2 2 3 3 1])',[2,size(dat.(zone).tri,1)*3])';
      % sort list
      edge = sortrows(sort(edge,2));
      % Mark outside edges
      mk = [1; sum(abs(diff(edge)),2)];
      mk = mk(2:end).*mk(1:end-1);
      % Store outside edges
      dat.(zone).edge = edge(find(mk~=0),:)';
    case('QUADRILATERAL'),
      dat.(zone).quad = reshape(fscanf(fid,'%f',dat.(zone).E*4),[4,dat.(zone).E])';
    otherwise
      error([mfilename ':help'],'unknown %s zone type %d',dat.(zone).name,dat.(zone).ET);
    end;

    % Plot
    if (plotem),
      switch(dat.(zone).ET),
      case('TRIANGLE'),
        HS(count) = trisurf(dat.(zone).tri,dat.(zone).X,dat.(zone).Y,dat.(zone).Z,dat.(zone).(vname{IV}));
        set(HS(count),'tag',dat.(zone).name,'edgecolor','none');
        HL{count} = plot3(dat.(zone).X(dat.(zone).edge),dat.(zone).Y(dat.(zone).edge),dat.(zone).Z(dat.(zone).edge),'k-');
        % Mirror
        %plot3(dat.(zone).X(dat.(zone).edge),-dat.(zone).Y(dat.(zone).edge),dat.(zone).Z(dat.(zone).edge),'k-');
        %HC(count) = trisurf(dat.(zone).tri,dat.(zone).X,-dat.(zone).Y,dat.(zone).Z,dat.(zone).(vname{IV}));
        %set(HC(count),'tag',['copy_of_' dat.(zone).name],'ButtonDownFcn','GTPopupText([],''none'');','edgecolor','none');
        drawnow;
      case('QUADRILATERAL'),
        HS(count) = trisurf(dat.(zone).quad,dat.(zone).X,dat.(zone).Y,dat.(zone).Z,dat.(zone).(vname{IV}));
        set(HS(count),'tag',dat.(zone).name,'edgecolor','none');
      otherwise
        error([mfilename ':help'],'unknown %s zone type %d',dat.(zone).name,dat.(zone).ET);
      end;
    end;

  end;
  ll = fgetl(fid);
end;
% Figure cleanup
if (plotem),
  hold off;
end;
% Close source data file
fclose(fid);

% Save data
if (saveem),
  [PATHSTR,NAME,EXT,VERSN] = fileparts(fname);
  if (~quiet), disp(sprintf('Saving compiled data to %s.mat',[PATHSTR filesep NAME])); end;
  save([PATHSTR filesep NAME],'dat');
  % Delete source data
  if (deleteem),
    if (~quiet), disp(sprintf('Deleting file %s',fname)); end;
    delete(fname);
  end;
end;
