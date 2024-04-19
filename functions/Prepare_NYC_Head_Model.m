function [dipin, grad, headmodel] = Prepare_NYC_Head_Model(channel, model_channel, sa, Res)

% This function prepares the NYC head model to become compatible with the
% Fieldtrip functions. It selects a subset of dipoles and their corresponding
% leadfield matrices. There are a subset of points available in the NYC model
% which indicates that those points should be considered as a whole in order to
% downsample the model. In each subset, the euclidean mean is calculated and 
% the nearest point to the mean is selected as a substitute for all of 
% the points in the corresponding subset. The leadfield matrix elements 
% containing the conductance of the path between these selected points and 
% the selected channels would be stored seperately. Because all of these
% dipoles are inside the brain, dipin.inside equals to ones(dipoles).

% INPUT: 
% 1. channel        = The EEG channels
% 2. model_channel  = The renamed NYC head model channels
% 3. sa             = The NYC head model .mat file
% 4. Res            = Indicated resolution

% OUTPUT:
% 1. dipin          = A struct containing the leadfields, the position of 
%                     the dipoles, the "inside" vector, and the dipole
%                     moments
% 2. grad,headmodel = [] --> These variables are set to [], because when
% the leadfield matrices are provided, there is no need to provide the
% headmodel in Fieldtrip.


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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find channel indices available in the NYC model
index_selected_channels = [];
for i = 1:length(channel)
    for j = 1:length(model_channel)
        if (strcmp(channel{i},model_channel{j}))
            index_selected_channels(i) = j;
            break;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the vertices of the selected resolution
Points_number = max(sa.(Res).tri(:,1));
Point_set = cell(Points_number,2); 
Center = cell(Points_number,1);
Selected_resolution_vc = zeros(Points_number,3);
Nearest = zeros(Points_number,1);

if (strcmp('cortex75K',Res) == 1)
    sa.(Res).in_to_cortex75K_eucl = 1:Points_number;
end

for i = 1:Points_number
    Point_set{i,1} = sa.cortex75K.vc((sa.(Res).in_to_cortex75K_eucl == i),:);
    Center{i,1} = mean(Point_set{i,1},1);
    Nearest(i) = dsearchn(Point_set{i,1},Center{i,1});
    Selected_resolution_vc(i,:) = Point_set{i,1}(Nearest(i),:);
    index_selected_points(i) = find((sa.cortex75K.vc(:,1) == Selected_resolution_vc(i,1))&(sa.cortex75K.vc(:,3) == Selected_resolution_vc(i,3))&(sa.cortex75K.vc(:,2) == Selected_resolution_vc(i,2)));
    dipin.leadfield{1,i} = squeeze(sa.cortex75K.V_fem(index_selected_channels,index_selected_points(i),:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining head model variables.
dipin.pos = Selected_resolution_vc;
dipin.mom = transpose(sa.cortex75K.normals(index_selected_points,:));
% All dipoles are inside the brain.
dipin.inside = ones(Points_number,1);
% By defining dipin.inside and dipin.leadfield, headmodel is not needed.
% By defining dipin.leadfield, grad is not used.
grad = [];
headmodel = [];

end