function [allowed, blocked] = civ_allowed_injection_intervals(zmin,zmax,zDoublet,zSinglet,velocityKms)
%CIV_ALLOWED_INJECTION_INTERVALS Exclude overlaps with original fitted lines.
% Each injected doublet member must be farther than velocityKms from every
% original doublet member and every original singlet 1548 feature. Intervals
% are continuous in redshift; no spectral pixels are modified or masked.
assert(isscalar(zmin)&&isscalar(zmax)&&isfinite(zmin)&&isfinite(zmax)&&zmin>=0&&zmax>zmin, ...
    'Invalid injection search window.');
assert(isscalar(velocityKms)&&isfinite(velocityKms)&&velocityKms>0&&velocityKms<299792.458, ...
    'Invalid exclusion velocity.');
assert(all(isfinite(zDoublet(:))&zDoublet(:)>=0) && all(isfinite(zSinglet(:))&zSinglet(:)>=0), ...
    'Invalid original absorber redshift.');
lam=[1548.1949462890625 1550.77001953125]; c=299792.458;
observed=[lam(1)*(1+zDoublet(:));lam(2)*(1+zDoublet(:));lam(1)*(1+zSinglet(:))];
blocked=zeros(numel(observed)*2,2); n=0;
for i=1:numel(observed)
    for j=1:2
        lo=max(zmin,observed(i)*(1-velocityKms/c)/lam(j)-1);
        hi=min(zmax,observed(i)*(1+velocityKms/c)/lam(j)-1);
        if hi>lo, n=n+1; blocked(n,:)=[lo hi]; end
    end
end
blocked=blocked(1:n,:);
if isempty(blocked), allowed=[zmin zmax]; return; end
blocked=sortrows(blocked,1); merged=blocked(1,:);
for i=2:size(blocked,1)
    if blocked(i,1)<=merged(end,2)
        merged(end,2)=max(merged(end,2),blocked(i,2));
    else
        merged(end+1,:)=blocked(i,:); %#ok<AGROW>
    end
end
blocked=merged; allowed=zeros(size(blocked,1)+1,2); n=0; left=zmin;
for i=1:size(blocked,1)
    if blocked(i,1)>left, n=n+1; allowed(n,:)=[left blocked(i,1)]; end
    left=max(left,blocked(i,2));
end
if left<zmax, n=n+1; allowed(n,:)=[left zmax]; end
allowed=allowed(1:n,:);
end
