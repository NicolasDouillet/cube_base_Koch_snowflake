function [V, T] = cubic_base_Koch_snowflake(nb_it, printable_ready, option_display)
%% cubic_base_Koch_snowflake : function to compute, display, and save
% the cubic base 3D Koch snowflake at any iteration / depth level.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%
% Syntax
%
% cubic_base_Koch_snowflake(nb_it);
% cubic_base_Koch_snowflake(nb_it, printable_ready);
% cubic_base_Koch_snowflake(nb_it, printable_ready, option_display);
% [V, T] = cubic_base_Koch_snowflake(nb_it, printable_ready, option_display);
%
%
% Description
%
% cubic_base_Koch_snowflake(nb_it) computes and displays the nb_it
% cubic base 3D Koch snowflake included in the unit sphere.
%
% cubic_base_Koch_snowflake(nb_it, printable_ready) prevents from
% creating non manifold edges when printable_ready is set to true /
% logical 1, and remove duplicated vertices and faces when it is set to
% *false / logical 0. In this latter case, the model is lighter (less
% vertices, less faces), but at the cost of non manifoldness.
%
% cubic_base_Koch_snowflake(nb_it, printable_ready, option_display)
% displays it when option_display is set to logical *true/1 (default),
% and doesn't when it is set to  logical false/0.
%
% [V,T] = cubic_base_Koch_snowflake(nb_it, printable_ready, option_display) saves
% the resulting vertex coordinates in the array V, and the triangulation in the array T.
%
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - printable_ready : either logical, true/*false or numeric 1/*0.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1
%
% Computes and displays the cubic base 3D Koch snowflake
% at iteration 3, with minimum vertex and face numbers
%
% cubic_base_Koch_snowflake(3);
%
%
% Example #2
%
% Computes and saves the cubic base 3D printable ready
% 3D Koch snowflake at iteration 3
%
% [V,T] = cubic_base_Koch_snowflake(4,true,false);


%% Input parsing
assert(nargin < 4,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    printable_ready = false;
    option_display = true;
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        assert(islogical(printable_ready) || isnumeric(printable_ready),'printable_ready parameter type must be either logical or numeric.');
        if nargin > 2
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
        else
            option_display = true;
        end
    else
        printable_ready = false;
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 3
    warning('%s triangles to display ! Make sure your graphic card has enough memory.',num2str(12*20^nb_it))    
end
warning('off');


%% Body
% Summits of original cube (living in the unit sphere R(O,1))
a = sqrt(3)/3;

V1 = [a a a];
V2 = [-a a a];
V3 = [-a -a a];
V4 = [a -a a];
V5 = -V1;
V6 = -V2;
V7 = -V3;
V8 = -V4;

C = cube(V1, V2, V3, V4, V5, V6, V7, V8);

% Loop on nb_it
p = 0;

while p ~= nb_it 
    
    new_C_array = repmat(C, [1 1 18]);
    
    for j = 1 : size(C,3)
        
        C_current = C(:,:,j);        
        [V_new, F_new] = split_cube(C_current);
        
        for m = 1:size(F_new,1)/6 %  18 = 108/6
            
            new_cube = cube(V_new(F_new(6*(m-1)+1,1),:),...
                            V_new(F_new(6*(m-1)+1,2),:),...
                            V_new(F_new(6*(m-1)+1,3),:),...
                            V_new(F_new(6*(m-1)+1,4),:),...
                            V_new(F_new(6*(m-1)+2,1),:),...
                            V_new(F_new(6*(m-1)+2,2),:),...
                            V_new(F_new(6*(m-1)+2,3),:),...
                            V_new(F_new(6*(m-1)+2,4),:)); % first two lines vertices in this order
                        
            new_C_array(:,:,18*(j-1) + m) = new_cube;
            
        end
        
    end
    
    C = new_C_array;        
    p = p+1;
    
end

% Squares to triangles conversion
[V,T] = squares2triangles(C);

if ~printable_ready
    
    % Remove duplicated vertices
    [V,T] = remove_duplicated_vertices(V,T);
    
    % Remove duplicated triangles
    T = unique(sort(T,2),'rows','stable');
    
end

% Display
if option_display
    
    cmap = [0 1 1];
    disp_Koch_snowflake(V,T,cmap);
    
end

end % cubic_base_Koch_snowflake


%% Cube structure computation subfunction
function [C] = cube(V1, V2, V3, V4, V5, V6, V7, V8)
%
% V1, V2, V3, V4, V5, V6, V7, V8 : line vectors of cube eight vertices coordinates
%
% Vn = [Vxn Vyn Vzn]

F1 = [1 2 3 4];
F2 = [5 6 7 8];
F3 = [1 4 8 5];
F4 = [2 1 5 6];
F5 = [2 3 7 6];
F6 = [3 4 8 7];

C = struct('vertex', [V1; V2; V3; V4; V5; V6; V7; V8], ...
           'facet', [F1; F2; F3; F4; F5; F6]);
       
end % cube


%% Split cube subfunction
function [V_new, F_new] = split_cube(C)
%
% Input
%
% C : cube structure
%
%
% Outputs
%
% V_new : The 56 newly created vertices -coordinates-
% F_new : 6 x 27 - 24 index line triplets of the 6 x 27 - 24 newly created facets


one_third = 1/3;
two_third = 2/3;

% Four possible values for each coordinate, X, Y, or Z
val1x = min(C.vertex(:,1)); % min value
val4x = max(C.vertex(:,1)); % max value
val2x = val1x + one_third * (val4x - val1x);
val3x = val1x + two_third * (val4x - val1x);

val1y = min(C.vertex(:,2)); % min value
val4y = max(C.vertex(:,2)); % max value
val2y = val1y + one_third * (val4y - val1y);
val3y = val1y + two_third * (val4y - val1y);

val1z = min(C.vertex(:,3)); % min value
val4z = max(C.vertex(:,3)); % max value
val2z = val1z + one_third * (val4z - val1z);
val3z = val1z + two_third * (val4z - val1z);


% Compute 56 new vertices coordinates | question de repérage...
% Regression order X, Y, Z for vertices

V_new = [val3x val4y val4z; % Top face
         val2x val4y val4z;
         val4x val3y val4z;
         val3x val3y val4z;
         val2x val3y val4z;
         val1x val3y val4z;
         val4x val2y val4z;
         val3x val2y val4z;
         val2x val2y val4z;
         val1x val2y val4z;
         val3x val1y val4z;
         val2x val1y val4z;...
         
         val4x val4y val3z; % Bottom first layer
         val3x val4y val3z; 
         val2x val4y val3z;
         val1x val4y val3z;
         val4x val3y val3z;
         val3x val3y val3z;
         val2x val3y val3z;
         val1x val3y val3z;
         val4x val2y val3z;
         val3x val2y val3z;
         val2x val2y val3z;
         val1x val2y val3z;
         val4x val1y val3z;
         val3x val1y val3z;
         val2x val1y val3z;
         val1x val1y val3z;...
         
         val4x val4y val2z; % Bottom second layer
         val3x val4y val2z; 
         val2x val4y val2z;
         val1x val4y val2z;
         val4x val3y val2z;
         val3x val3y val2z;
         val2x val3y val2z;
         val1x val3y val2z;
         val4x val2y val2z;
         val3x val2y val2z;
         val2x val2y val2z;
         val1x val2y val2z;
         val4x val1y val2z;
         val3x val1y val2z;
         val2x val1y val2z;
         val1x val1y val2z;...
         
         val3x val4y val1z; % Bottom face
         val2x val4y val1z;
         val4x val3y val1z;
         val3x val3y val1z;
         val2x val3y val1z;
         val1x val3y val1z;
         val4x val2y val1z;
         val3x val2y val1z;
         val2x val2y val1z;
         val1x val2y val1z;
         val3x val1y val1z;
         val2x val1y val1z]; 

% 6 x 27 - 24 = 138 new facets but actually 6 x 18 = 108 new facets
% /_!_\ Counter clockwise sorted for squares /_!_\

% General model  for a (a b c d e f g h) cube
% a b c d
% e f g h
% a d h e
% b a e f
% c b f g
% d c g h

F_new = [1 2 5 4; % Top layer top cross cube
         14 15 19 18;
         1 4 18 14;
         2 1 14 15;
         5 2 15 19;
         4 5 19 18;...
         
         3 4 8 7; % Top layer right cross cube
         17 18 22 21;
         3 7 21 17;
         4 3 17 18;
         8 4 18 22;
         7 8 22 21;...
         
         4 5 9 8; % Top layer centre cross cube
         18 19 23 22;
         4 8 22 18;
         5 4 18 19;
         9 5 19 23;
         8 9 23 22;...
         
         5 6 10 9; % Top layer left cross cube
         19 20 24 23;
         5 9 23 19;
         6 5 19 20;
         10 6 20 24;
         9 10 24 23;...
         
         8 9 12 11; % Top layer bottom cross cube
         22 23 27 26;
         8 11 26 22;
         9 8 22 23;
         12 9 23 27;
         11 12 27 26;...         
        
         13 14 18 17; % Middle layer top right square cube
         29 30 34 33;
         13 17 33 29;
         14 13 29 30;
         18 14 30 34;
         17 18 34 33;...         
         
         14 15 19 18; % Middle layer top cross cube
         30 31 35 34;
         14 18 34 30;
         15 14 30 31;
         19 15 31 35;
         18 19 35 34;...
         
         15 16 20 19; % Middle layer top left square cube
         31 32 36 35;
         15 19 35 31;
         16 15 31 32;
         20 16 32 36;
         19 20 36 35;...                          
         
         17 18 22 21; % Middle layer right cross cube
         33 34 38 37;
         17 21 37 33;
         18 17 33 34;
         22 18 34 38;
         21 22 38 37;...
         
         19 20 24 23; % Middle layer left cross cube
         35 36 40 39;
         19 23 39 35;
         20 19 35 36;
         24 20 36 40;
         23 24 40 39;...
         
         21 22 26 25; % Middle layer bottom right square cube
         37 38 42 41; 
         21 25 41 37;
         22 21 37 38;
         26 22 38 42;
         25 26 42 41;...
          
         22 23 27 26; % Middle layer bottom cross cube
         38 39 43 42;
         22 26 42 38;
         23 22 38 39;
         27 23 39 43;
         26 27 43 42;...
          
         23 24 28 27; % Middle layer bottom left square cube
         39 40 44 43;
         23 27 43 39;
         24 23 39 40;
         28 24 40 44;
         27 28 44 43;...
          
         
         30 31 35 34; % Bottom layer top cross cube
         45 46 49 48;
         30 34 48 45;
         31 30 45 46;
         35 31 46 49;
         34 35 49 48;...
         
         33 34 38 37; % Bottom layer right cross cube
         47 48 52 51;
         33 37 51 47;
         34 33 47 48;
         38 34 48 52;
         37 38 52 51;...
         
         34 35 39 38; % Bottom layer centre cross cube
         48 49 53 52;
         34 38 52 48;
         35 34 48 49;
         39 35 49 53;
         38 39 53 52;...         
         
         35 36 40 39; % Bottom layer left cross cube
         49 50 54 53;
         35 39 53 49;
         36 35 49 50;
         40 36 50 54;
         39 40 54 53;...         
         
         38 39 43 42; % Bottm layer bottom cross cube
         52 53 56 55;
         38 42 55 52;
         39 38 52 53;
         43 39 53 56;
         42 43 56 55];
     
end % split_cube


%% Squares to triangles conversion subfunction
function [V, T] = squares2triangles(C)
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.
%
% Split struct array into two arrays : vertices & facets


S = size(C,3);
V = zeros(8*S,3);
T = zeros(12*S,3);

for k = 1:S
    
    for i = 1:size(C(:,:,k).vertex,1)
        
        V(8*(k-1)+i,:) = C(:,:,k).vertex(i,:);
        
    end
    
    for j = 1:size(C(:,:,k).facet,1) % 6
        
        a = C(:,:,k).facet(j,1) + 8*(k-1);
        b = C(:,:,k).facet(j,2) + 8*(k-1);
        c = C(:,:,k).facet(j,3) + 8*(k-1);
        d = C(:,:,k).facet(j,4) + 8*(k-1);
        
        T1 = sort([a b c]);
        T2 = sort([a d c]);
        
        T(12*(k-1)+2*(j-1)+1,:) = T1;
        T(12*(k-1)+2*j,:) = T2;
        
    end
    
end

% Remove duplicated triangles
T = unique(T,'rows','stable');

end % squares2triangles


%% Display subfunction
function [] = disp_Koch_snowflake(V, T, cmap)
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.

figure;
set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
trisurf(T,V(:,1),V(:,2),V(:,3),'EdgeColor',cmap), shading interp, hold on;
colormap(cmap);
axis square, axis equal, axis tight, axis off;
grid off;
ax = gca;
ax.Clipping = 'off';
camlight left;
view(-45,35);

end % disp_Koch_snowflake


%% Remove duplicated vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices