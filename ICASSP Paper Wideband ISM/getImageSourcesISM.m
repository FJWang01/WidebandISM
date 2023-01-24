function Image = getImageSourcesISM(Reflection_order, Room, Source, Receiver, fs, N, c, time_vec, LPF_order, LPF_fn)  
% This function calculates the properties of the image sources including
% the original source
%
% Input
% Reflection_order - reflection order, scalar
%
% Room - structure with fields
%        Room.size - room dimension, 1 by 3 vector
%        Room.b - wall reflection coefficients, 1 by 6 vector
%
% Source - structure with fields
%          Source.GL_Os - location of the source w.r.t. the origin of the room
%          Source.GL_SH_coeffs - SH coefficients of the directivities of 
%                                the spherical wave fronts
%          Source.tau - time instances when spherical wave fronts are
%                       emitted, must be a row vector
%
% Receiver - structure with fields
%            Receiver.GL_Or - location of the receiver origin w.r.t. the
%                          origin of the room
%            Receiver.R - radius of the observation sphere
% fs - sampling frequency, in Hz
% N - spherical harmonic truncation order of the observed signal on the
%     observation sphere
% c - speed of sound, in metres per second
% time_vec - time sample vector, in seconds
% LPF_order - order of the low-pass filter
% LPF_fn - cut-off frequency of the low-pass filter, in Hz
%
% Output
% Image - structure with fields
%         %% GL stands for global
%         Image.GL_x - x coordinates of the image sources w.r.t.
%                                 the origin of the room
%         Image.GL_y - y coordinates of the image sources w.r.t.
%                                 the origin of the room
%         Image.GL_z - z coordinates of the image sources w.r.t.
%                                 the origin of the room
%
%         %% Coordinates of the image sources w.r.t. the receiver origin Or
%         Image.Or_x - x coordinates of the image sources w.r.t. the
%                      receiver origin
%         Image.Or_y - y coordinates of the image sources w.r.t. the
%                      receiver origin
%         Image.Or_z - z coordinates of the image sources w.r.t. the
%                      receiver origin
%
%         Image.Or_d - radii of the image sources w.r.t. the receiver
%                      origin
%         Image.Or_elev - elevations of the image sources w.r.t. the
%                         receiver origin
%         Image.Or_azim - azimuths of the image sources w.r.t. the receiver
%                         origin
%
%         Image.refl_coeff - reflection coefficients of the image sources
%
%         Image.scaling - scalings of the directivity patterns of the 
%                          image sources due to reflection (1 or -1)
%         Image.GL_degree - degrees of the directivities of the spherical
%                            wave fronts emitted by the image sources 
%         Image.GL_order - orders of the directivities of the spherical
%                           wave fronts emitted by the image sources 
%         Image.num_of_images - number of image sources for the given 
%                               reflection order

%%  Input Checking  
% Reflection_order
validateattributes(Reflection_order, {'double'}, {'integer','scalar'});
                    
% Room
if ~isstruct(Room)
    error('Expected input Room to be a struct');
else
    validateattributes(Room.size,{'double'}, {'size',[1,3]});
    validateattributes(Room.b,{'double'}, {'size',[1,6],'>=',0,'<=',1});
end

% Source
if ~isstruct(Source)
    error('Expected input Source to be a struct');
else    
    validateattributes(Source.GL_Os,{'double'}, {'size',[1,3],'positive'}); 
    validateattributes(Source.tau, {'double'}, {'row', 'nonnegative'});
end

if ~isequal(size(Source.GL_SH_coeffs, 2), numel(Source.tau))
    error(['The number of spherical wave fronts must match the number of ' ...
        'time instances when they are emitted']);
else
    % do nothing
end

% Receiver
if ~isstruct(Receiver)
    error('Expected input Receiver to be a struct');
else    
    validateattributes(Receiver.GL_Or,{'double'}, {'size',[1,3],'positive'});  
    validateattributes(Receiver.R,{'double'}, {'scalar'}); 
end

%% Generate image source 
% 1. Triple indexing (p1, p2, p3) [1 by 8]
%    p1, p2, p3 can be either 0 or 1

%   If p1 = 1, then an image source in the x-axis direction is active
%   If p2 = 1, then an image source in the y-axis direction is active
%   If p3 = 1, then an image source in the z-axis direction is active
p1 = [0, 1, 0, 1, 0, 1, 0, 1];    
p2 = [0, 0, 1, 1, 0, 0, 1, 1];
p3 = [0, 0, 0, 0, 1, 1, 1, 1];

% 2. m_x, m_y, m_z indexing
%    Reflection_order is the reflection order
%    Calculate the maximum of m_x or m_y or m_z
max_val_of_m = ceil((Reflection_order + 1)/2);

%    Valid values of m_x, m_y, m_z
valid_val_of_m = -max_val_of_m:max_val_of_m;

%    All valid combinations from valid values of m_x, m_y, m_z
[comb_of_m, ~] = permn(valid_val_of_m, 3);

%   Allocate m_x, m_y, m_z
m_x = comb_of_m(:, 1);
m_y = comb_of_m(:, 2);
m_z = comb_of_m(:, 3);

% change column vector to row vector
m_x = m_x.'; 
m_y = m_y.';
m_z = m_z.';

% 3. Generate all combinations of m_x, m_y, m_z, p1, p2, p3
%    Repeat the whole m_x, m_y, m_z vectos 8 times
m_x_rep = repmat(m_x, 1, numel(p1));
m_y_rep = repmat(m_y, 1, numel(p1));
m_z_rep = repmat(m_z, 1, numel(p1));

%    Repeat each element of p1, p2, p3 by numel(m_x) times
p1_rep = repelem(p1, numel(m_x) * ones(1, numel(p1)));
p2_rep = repelem(p2, numel(m_x) * ones(1, numel(p2)));
p3_rep = repelem(p3, numel(m_x) * ones(1, numel(p3)));

% 4. Find  m_x, m_y, m_z, p1, p2, p3 combinations that are less than or
%    equal to Reflection_order
comb_order = abs(2*m_x_rep - p1_rep) + abs(2*m_y_rep - p2_rep) + abs(2*m_z_rep - p3_rep);
valid_images = comb_order <= Reflection_order;

% 5. Calculate the locations of image sources including the original source
Image.GL_x = (1-2*p1_rep(valid_images))*Source.GL_Os(1) + 2*m_x_rep(valid_images)*Room.size(1);
Image.GL_y = (1-2*p2_rep(valid_images))*Source.GL_Os(2) + 2*m_y_rep(valid_images)*Room.size(2);
Image.GL_z = (1-2*p3_rep(valid_images))*Source.GL_Os(3) + 2*m_z_rep(valid_images)*Room.size(3);

%    Count the number of image sources including the original source
Image.num_of_images = numel(Image.GL_x);

% 6. Find the vectors pointing from the receiver origin to the image sources
%    Cartesian 
Image.Or_x = Image.GL_x - Receiver.GL_Or(1);
Image.Or_y = Image.GL_y - Receiver.GL_Or(2);
Image.Or_z = Image.GL_z - Receiver.GL_Or(3);

%    Spherical 
Image.Or_d = sqrt(Image.Or_x.^2 + Image.Or_y.^2 + Image.Or_z.^2);
Image.Or_elev = atan2(sqrt(Image.Or_x.^2 + Image.Or_y.^2), Image.Or_z); 
Image.Or_azim = wrapTo2Pi(atan2(Image.Or_y, Image.Or_x));

% 7. Find the reflection coefficient of each image source
Image.refl_coeff = (  Room.b(1).^abs(m_x_rep(valid_images) - p1_rep(valid_images)) ...
                    .* Room.b(2).^abs(m_x_rep(valid_images))  ...
                    .* Room.b(3).^abs(m_y_rep(valid_images) - p2_rep(valid_images)) ...
                    .* Room.b(4).^abs(m_y_rep(valid_images)) ...
                    .* Room.b(5).^abs(m_z_rep(valid_images) - p3_rep(valid_images)) ...
                    .* Room.b(6).^abs(m_z_rep(valid_images)) ...
                    );
                
% 8. Find the degree and order of each image source
% Image.scalings = (-1).^((p2_rep(valid_images) + p3_rep(valid_images))*Source.order...
%                     + p3_rep(valid_images)*Source.degree);
% Image.degrees = Source.degree*ones(1, Image.num_of_image_sources);
% Image.orders = (-1).^(p1_rep(valid_images) + p2_rep(valid_images)) * Source.order;

% max SH order of the spherical wave front emitted by the source
max_GL_order = sqrt(size(Source.GL_SH_coeffs, 1))-1;
[n_all, m_all, ~] = getDegreeOrderPairs(max_GL_order);
Image.scalings = (-1).^(kron((p2_rep(valid_images) + p3_rep(valid_images)),m_all)...
                    + kron(p3_rep(valid_images), n_all));
Image.degrees = repmat(n_all, 1, Image.num_of_images);
Image.orders = kron((-1).^(p1_rep(valid_images) + p2_rep(valid_images)) , m_all);

%% Determine if there is the need to flip opposite order
order_scaling = (-1).^(p1_rep(valid_images) + p2_rep(valid_images));
% Pick -1 where -1 means SH coefficients of the opposite order need to be
% swapped
swap_mat = order_scaling < 0;

%% Apply scaling due to reflections and flip opposite order when desired
% size(Image.GL_SH_coeff) = [size(Source.GL_SH_coeffs, 1), numel(valid_images), size(Source.GL_SH_coeffs, 2)]
% image sources are the 2nd dimension
% tau is the 3rd dimension
Image.GL_SH_coeff = image_source_swap_orders(Source.GL_SH_coeffs, swap_mat, Image.scalings);
Image.GL_SH_coeff = Image.GL_SH_coeff .* Image.refl_coeff;
%% Rotate directivities of spherical wave fronts to form axis symmetric case
azim_inv_rot = wrapTo2Pi(-Image.Or_azim - pi);
Wigner_D_cell = arrayfun(@(x, y) WignerDM(max_GL_order, pi, x, y), Image.Or_elev, azim_inv_rot, 'UniformOutput', false);
% size(Wigner_D_mat) = [(max_GL_order+1)^2, size(Source.GL_SH_coeffs, 1), numel(valid_images)]
Wigner_D_mat = cat(3, Wigner_D_cell{:});

% size(im_GL_SH_coeff_shu) = [size(Source.GL_SH_coeffs, 1), size(Source.GL_SH_coeffs, 2), numel(valid_images)]
% tau is the 2nd dimension
% image sources are the 3rd dimension
im_GL_SH_coeff_shu = permute(Image.GL_SH_coeff, [1, 3, 2]);

% size(im_sym_SH_coeff_shu) = [size(Source.GL_SH_coeffs, 1), size(Source.GL_SH_coeffs, 2), numel(valid_images)]
im_sym_SH_coeff_shu = pagemtimes(Wigner_D_mat, im_GL_SH_coeff_shu);

% size(im_sym_SH_coeff) = [size(Source.GL_SH_coeffs, 1), numel(valid_images), size(Source.GL_SH_coeffs, 2)]
% image sources are the 2nd dimension
% tau is the 3rd dimension
im_sym_SH_coeff = permute(im_sym_SH_coeff_shu, [1, 3, 2]);
Image.sym_SH_coeff_cell = num2cell(im_sym_SH_coeff);
%% Calculate the SH coefficients of the observed signal 
% SH coefficients are for the axis symmetric case (post rotation needed)
% Rearrange variable dimension to make sure:
% n, m are the 1st dimension
% image sources are the 2nd dimension
% tau is the 3rd dimension
Or_d_mat = repmat(Image.Or_d, numel(n_all), 1, numel(Source.tau));
n_all_mat = repmat(n_all, 1, numel(Image.Or_d), numel(Source.tau));
m_all_mat = repmat(m_all, 1, numel(Image.Or_d), numel(Source.tau));
tau_mat = reshape(Source.tau, 1, 1, numel(Source.tau));
tau_mat = repmat(tau_mat, numel(n_all), numel(Image.Or_d));

% SH_coeff = higher_order_source_on_z_axis_delayed(r, r_s, v, u, N, time_vec, c, tau)
im_Obv_sym_SH_coeff_ind = arrayfun(@(x, y, z, p) higher_order_source_on_z_axis_delayed(Receiver.R, x, y, z, N, time_vec, c, p),...
                                Or_d_mat, n_all_mat, m_all_mat, tau_mat, 'UniformOutput', false);

Image.Obv_sym_SH_coeff = cellfun(@(x, y) x.*y, Image.sym_SH_coeff_cell, im_Obv_sym_SH_coeff_ind, 'UniformOutput', false);
%% Post rotation
Wigner_D_cell_po = arrayfun(@(x, y) WignerDM(N, x, y, 0), Image.Or_azim, Image.Or_elev, 'UniformOutput', false); 
Wigner_D_cell_po_rep = repmat(Wigner_D_cell_po, numel(n_all), 1, numel(Source.tau));

Image.Obv_GL_SH_coeff = cellfun(@(x, y) x*y, Wigner_D_cell_po_rep, Image.Obv_sym_SH_coeff, 'UniformOutput', false);

% Sum contribution of all image sources
Image_cat = cat(3, Image.Obv_GL_SH_coeff{:});
Image.Obv_GL_SH_coeff_sum = sum(Image_cat, 3);

%% Design FIR filter
Wn = LPF_fn/(fs/2); % normalised cut-off frequency
b = fir1(LPF_order, Wn);

% visualise filter frequency response
% h = fvtool(b, 1);
% get(h) % display all fields of h
% h.Fs = fs;
% h.FrequencyRange = '[0, Fs/2)';

SH_coeff_scaled_and_smoothed = filter(b, 1, Image.Obv_GL_SH_coeff_sum, [], 2);
Image.Obv_GL_SH_coeff_final = circshift(SH_coeff_scaled_and_smoothed, -LPF_order/2, 2);
end % function
