%% cube_base_Koch_snowflake
%
% Function to compute, display, and save the cube base 
% 3D Koch snowflake at any iteration.
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
%
%% Syntax
%
% cube_base_Koch_snowflake;
% 
% cube_base_Koch_snowflake(nb_it);
%
% cube_base_Koch_snowflake(nb_it, printable_ready);
%
% cube_base_Koch_snowflake(nb_it, printable_ready, option_display);
%
% [V, T] = cube_base_Koch_snowflake(nb_it, printable_ready, option_display);
%
%% Description
%
% cube_base_Koch_snowflake computes and displays the 3-
% cube base 3D Koch snowflake included in the unit sphere.
%
% cube_base_Koch_snowflake(nb_it) computes and displays the nb_it
% cube base 3D Koch snowflake included in the unit sphere.
%
% cube_base_Koch_snowflake(nb_it, printable_ready) prevents from
% creating non manifold edges when printable_ready is set to true /
% logical 1, and remove duplicated vertices and faces when it is set to
% *false / logical 0. In this latter case, the model is lighter (less
% vertices, less faces), but at the cost of non manifoldness.
%
% cube_base_Koch_snowflake(nb_it, printable_ready, option_display)
% displays it when option_display is set to logical *true/1 (default),
% and doesn't when it is set to  logical false/0.
%
% [V,T] = cube_base_Koch_snowflake(nb_it, printable_ready, option_display) saves
% the resulting vertex coordinates in the array V, and the triangulation in the array T.
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/73247-n-level-3d-koch-snowflake Koch_snowflake> |
%   <https://fr.mathworks.com/matlabcentral/fileexchange/73447-n-level-3d-sierpinski-menger-sponge Sierpinski_Menger_sponge> |
%
%% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - printable_ready : either logical, true/*false or numeric 1/*0.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the cube base 3D Koch snowflake at iteration 3, with minimum vertex and face numbers

cube_base_Koch_snowflake(3);

%% Example #2
% Computes and saves the cube base 3D printable ready 3D Koch snowflake at iteration 3

[V,T] = cube_base_Koch_snowflake(4,true,false);