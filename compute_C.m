function ret = compute_C(phi)
%M * L * F, col major
global LAMBDAS D F SHAPE_M M


phi = reshape(phi, [1, 1, numel(phi)]);
C_M1 = reshape(0:SHAPE_M(1)-1, [SHAPE_M(1), 1]);
C_M2 = reshape(0:SHAPE_M(2)-1, [1, SHAPE_M(2)]);

C_ML = bsxfun(@plus, bsxfun(@times, C_M1, sin(phi)), bsxfun(@times, C_M2, cos(phi)));
C_LAMBDA = reshape(LAMBDAS, [1,1,1,F]);
ret = exp(1j*2*pi*D.*bsxfun(@rdivide, C_ML, C_LAMBDA));
% compress M for one dimension
ret = reshape(ret, [M, numel(phi), F]);
end