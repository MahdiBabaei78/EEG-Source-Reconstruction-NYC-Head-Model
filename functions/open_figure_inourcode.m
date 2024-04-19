function h = open_figure_inourcode(cfg)

% Original file: open_figure (Fieldtrip toolbox)
% Modified file: open_figure_inourcode
% Date: April 15, 2024
% Copyright (C) 2024 Mahdi Babaei

% This file has been modified from its original version. The following 
% changes have been made:
% 1. To create a UI figure, I modified the Fieldtrip code to create a figure
% at a specified position using the uifigure MATLAB function. (Line 66)

% Copyright (C) 2024-, Mahdi Babaei and Bjorn Erik Juel
%
% This file is part of EEG Source Reconstruction - NYC head model, see 
% https://github.com/MahdiBabaei78/EEG-Source-Reconstruction-NYC-Head-Model
% for the documentation and details.
%
%    EEG-Source-Reconstruction-NYC-Head-Model is a free software: 
%    you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EEG-Source-Reconstruction-NYC-Head-Model is distributed in the hope
%    that it will be useful, but WITHOUT ANY WARRANTY; without even the
%    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
%    PURPOSE.  See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EEG-Source-Reconstruction-NYC-Head-Model. If not, 
%    see <http://www.gnu.org/licenses/>.
%
% $Id$



% OPEN_FIGURE is a helper function to open a figure with some specific settings
% consistent over all FieldTrip functions that do plotting and/or that show a
% graphical user interface.
%
% See also GCA, GCF, GROOT, 

cfg.figure      = ft_getopt(cfg, 'figure', 'yes');
cfg.visible     = ft_getopt(cfg, 'visible', 'yes');
cfg.position    = ft_getopt(cfg, 'position'); % use the default position
cfg.renderer    = ft_getopt(cfg, 'renderer'); % let MATLAB decide on the default
cfg.figurename  = ft_getopt(cfg, 'figurename');
cfg.title       = ft_getopt(cfg, 'title');

switch cfg.figure
  case {'new', 'yes'}
    figopt = {};
    if ~istrue(cfg.visible)
      figopt = ft_setopt(figopt, 'visible', 'off');
    end
    if ~isempty(cfg.position)
      figopt = ft_setopt(figopt, 'position', cfg.position);
    end
    
    % check whether there is already a figure open
    h = get(0, 'CurrentFigure');

    if isempty(h)
      % there is no figure open yet, make a new one
      %%%%%%%%%%%%%%%%%%%% Modification %%%%%%%%%%%%%%%%%%
      h = uifigure("Position",[50 50 800 700]);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isempty(get(h, 'Children'))
      % there is an empty figure open, use that and update it
      if ~isempty(figopt)
        set(h, figopt{:});
      end
    else
      % there is another figure open, make a new one
      h = figure(figopt{:});
    end
    
  case {'clf', 'no'}
    % use the current figure and clear it
    % this will open a new one if there is no figure yet
    h = clf;

  case 'gca'
    % use the current axes but do not clear them
    % this will open a new figure if there is no figure yet
    h = gca;

  case 'gcf'
    % use the current figure but do not clear it
    % this will open a new figure if there is no figure yet
    h = gcf;
  
  otherwise
    % assume that it specifies a figure handle
    h = cfg.figure;
end

assert(ishandle(h), 'failed to open figure');

if ~isempty(cfg.figurename)
  % this appears as the name of the window
  set(h, 'name', cfg.figurename);
end

if ~isempty(cfg.renderer)
  % set figure renderer
  set(h, 'Renderer', cfg.renderer);
end

if ~isempty(cfg.title)
  % this appears above the axes
  title(cfg.title);
end
