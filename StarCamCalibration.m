% Geoffrey Palo
% EN.525.665 Machine Perception
% Final Project - Self-calibration of Star Cameras
% Algorithms from "Autonomous Star Camera Calibration and Spacecraft 
% Attitude Determination", Madhumita Pal & M. Seetharama Bhat 2017

clear

% point correspondences (calculated centroids using other software)
load('p_init1.mat');

% ASTHROS camera previously calculated principal axis offset
u0_goal = 83.5658;
v0_goal = 55.5521;

s_pix = 4.5e-6; % Pixel size from sensor datasheet
f = 135e-3; % Lens focal length
alpha_goal = f/s_pix; % assumed alpha/beta values
beta_goal = alpha_goal;

% Point correspondences of stars in ASTHROS moon images
p(:,:,1) = p_init(1:4,:);
p(:,:,2) = p_init(5:8,:);
p(:,:,3) = p_init(9:12,:);
p(:,:,4) = p_init(13:16,:);

% Convert pixel map values to u-v axes (subtract midpoint)
p(:,1,:) = p(:,1,:) - 4096/2;
p(:,2,:) = 4096/2 - p(:,2,:);

% A*m = 0, where m is the vector of M elements (derivation in Pal & Bhat)
A = zeros(8,9,size(p,3)-1);

for idx = 2:size(p,3)
    A_temp = [[p(1,:,1) 1 zeros(1,3) -p(1,1,idx)*[p(1,:,1) 1]];...
        [zeros(1,3) p(1,:,1) 1 -p(1,2,idx)*[p(1,:,1) 1]]];
    % This computes A for every POINT CORRESPONDENCE in one IMAGE PAIR
    for jdx = 2:size(p,1)
        A_temp = [A_temp; [[p(jdx,:,1) 1 zeros(1,3) -p(jdx,1,idx)*[p(jdx,:,1) 1]];...
            [zeros(1,3) p(jdx,:,1) 1 -p(jdx,2,idx)*[p(jdx,:,1) 1]]]];
    end
    % This creates the array A which contains N-1 8x9 matrices (where N is
    % number of images)
    A(:,:,idx - 1) = A_temp;
end

% Vector of homography matrix M elements
m = zeros(9,size(p,3)-1);

% Singular value decomposition of AT*A, m is the eigenvector corresponding
% to the smallest eigenvalue
for idx = 1:size(A,3)
    [V,D] = eig(A(:,:,idx)'*A(:,:,idx));
    [U,S,Vs] = svd(A(:,:,idx)'*A(:,:,idx));
    m(:,idx) = Vs(:,end);
end

M = zeros(3,3,size(m,2));
% Reshape and normalize m to M
for idx = 1:size(m,2)
    M(:,:,idx) = [m(1:3,idx)'; m(4:6,idx)'; m(7:9,idx)'];
    M(:,:,idx) = 1/M(3,3,idx)*M(:,:,idx);
end

Q = [];

% Create matrix Q, to solve K*K'*M^(-1)' = M*K*K'
for kdx = 1:size(M,3)
    M_invT(:,:,kdx) = inv(M(:,:,kdx))';
    Q = [Q; kron(eye(3),M(:,:,kdx)) - kron(M_invT(:,:,kdx),eye(3))];
end

% Singular value decomposition of QT*Q, k is the eigenvector corresponding
% to the smallest eigenvalue
[U,S,V] = svd(Q'*Q);
k = V(:,end);
KKt = reshape(k,[3,3]); 
KKt = 1/KKt(3,3) * KKt;

% Extract camera parameters from matrix K*K'
u0 = KKt(3,1);
v0 = KKt(3,2);
alpha = sqrt(KKt(1,1) - u0^2);
beta = sqrt(KKt(2,2) - v0^2);

% Starting to test k1 values

% Range of k1 (from Pal & Bhat)
r = sqrt(2048^2+2048^2);
k1Max = 2/r^2;

% For setting parameters to previously calculated, if desired
% alpha = alpha_goal;
% beta = beta_goal;
% u0 = u0_goal;
% v0 = v0_goal;

% Fifth point correspondence; found in images 1, 2, and 4
p1 = [2453.313477; 245.677963; 1];
p_actual(:,1) = [2405.280762; 350.052643; 1];
p_actual(:,2) = [2133.58374; 895.180969; 1];

% Convert pixel map values to u-v axes (subtract midpoint)
p1 = p1 - [4096/2; 4096/2; 0];
p_actual = [4096/2; 4096/2; 0] - p_actual;

for kdx = 1:size(p_actual,2)
    % Estimated point location using homography matrix
    if(kdx ~= 2)
        p_est = M(:,:,kdx)*p1;
    else
        p_est = M(:,:,kdx+1)*p1;
    end
    u = p_est(1);
    v = p_est(2);
    u_actual = p_actual(1,kdx);
    v_actual = p_actual(2,kdx);
    k1 = -k1Max:10^-5*k1Max:k1Max; % test values for k1
    for idx = 1:size(k1,2)
        % Calculate distorted point location with test k1 values
        ud = u + (u - u0)*k1(idx)*((u-u0)^2/alpha^2+(v-v0)^2/beta^2);
        vd = v + (v - u0)*k1(idx)*((u-u0)^2/alpha^2+(v-v0)^2/beta^2);
        % Distances of distorted (calculated) and actual points
        dist_d = sqrt(ud^2 + vd^2);
        dist_actual = sqrt(u_actual^2 + v_actual^2);
        % Squared error
        err(idx,kdx) = (dist_d - dist_actual)^2;
    end
end

[min_err, I] = min(err);

plot(k1,err(:,1),k1,err(:,2));
grid on;
xlabel('k1');
ylabel('Least square error');
legend('M1','M3');

alpha
beta
u0
v0
k1_M1 = k1(I(1))
k1_M2 = k1(I(2))