%     %% RUNNER UNTUK 1-POINT

clearvars
% tStart = tic;
load square_5nm

figure()
patch('Vertices',p,'Faces',t,'FaceColor','g','EdgeColor','k');
for m=1:length(t)
    N=t(m,1:3);
    X(:,m) = p(N,1)*10^9;
    Y(:,m) = p(N,2)*10^9;
    Z(:,m) = p(N,3)*10^9;
end
rotate3d
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')

Triangulation
SingularSpecial
SingularClose2_vect
% Surf_current_vectorized
farfield_optical_constants
% 
% field_farfield
%
