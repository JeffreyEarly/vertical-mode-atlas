function N2function = SmoothN2(N2,z,zLim)
K = 3;
constraints = [];

% These functions are use to transform to and from log-log space of (z,N2)
sz = @(z) flip(log10(1-z));
%zs = @(s) flip(1-10.^s);
gN2 = @(N2) flip(log10(N2));
N2g = @(g) flip(10.^g);

% Group the abyssal (<2500m) ocean measurements together
DF = 2;
abyssalDF = 5;
abyssalDepth = -2000;
abyssalIndex = find(z<abyssalDepth,1,'last');
lastIndex = 1;
if ~isempty(abyssalIndex) && (abyssalIndex-lastIndex)>DF
    % Makes sure there are a minimum of abyssalDF points between knot points
    if (abyssalIndex-lastIndex)/abyssalDF >= 2
        stride=floor((abyssalIndex-lastIndex)/floor((abyssalIndex-1)/abyssalDF));
        z_knot = [min(zLim); z(1+stride:stride:abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
    else
        z_knot = [min(zLim); z(abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
    end
    % Force the top and bottom to only allow for linear splines (no second derivative).
    if K > 2
        constraints = struct('t',sz([z_knot(1)+(z_knot(2)-z_knot(1))/2; z(end-1)+(z(end)-z(end-1))/2]),'D',[2;2]);
    end
else
    z_knot = [min(zLim); z(1+DF:DF:end-DF); max(zLim)];
    % Force the top to only allow for linear splines (no second derivative).
    if K > 2
        constraints = struct('t',sz(z(end-1)+(z(end)-z(end-1))/2),'D',2);
    end
end

% Now go ahead and do the correct boundary knot repeats/removals
z_knot = InterpolatingSpline.KnotPointsForPoints(z_knot,K,1);



% The error in log space is mostly proportional to the magnitude of the function.
% The factor of 4 includes the sqrt, but also just works well.
sigma = flip(10.^(abs(log10(N2/max(N2)))/4 ));


spline = SmoothingSpline(sz(z),gN2(N2),NormalDistribution(1),'t_knot',sz(z_knot),'K',K,'constraints',constraints,'sigma',sigma);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromGCV);
if spline.isConstrained == 1
    DiffMatrix = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(spline.T,spline.t,1);
    a = DiffMatrix*spline.x;
    a_rms = sqrt(mean(a.*a));
    spline.lambda = 1/(a_rms^2);
    spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromGCV);
end

if spline.isConstrained == 1
    spline.lambda = 0;
    fprintf('Spline still constrained! Setting smoothing parameter to zero.\n');
end


N2function = @(zz) N2g(spline(sz(zz)));
end