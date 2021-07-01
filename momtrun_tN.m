function [mv,var,nk]=momtrun_tN(m,s,a,b)
% function evaluates moments of truncated normal distribution
% keyboard
alf	= (a-m)./(sqrt(2)*s);
bet	= (b-m)./(sqrt(2)*s);

% Evaluation has sense for stable only
istab	= find((alf<=3) & (bet>=-3));
astab	= find(alf>3);
bstab	= find(bet<-3);
% keyboard
%%%%%%%% STABLE
if length(istab>0)
	alf	= alf(istab);
	bet	= bet(istab);

	pom	= (erf(bet) - erf(alf))* sqrt(pi/2);
	gam	= (exp(-bet.^2) - exp(-alf.^2)) ./ pom;
	del	= (b.*exp(-bet.^2) - a.*exp(-alf.^2)) ./ pom;

	mv(istab)	= m(istab) - s(istab).*gam;
	var(istab)	= s(istab).^2 + m(istab).*mv(istab)' - s(istab).*del;
	nk(istab)	= s(istab) .* pom;
end
%%%%%%%% UNSTABLE
if length(astab)>0
	ma	= m(astab);
	sa	= s(astab);
	mv(astab)	= a - sa.^2./ma - (a*sa.^2)./(ma.^2) - (a^2*sa.^2-2*sa.^4)./(ma.^3);
	var(astab)	= a^2 - (2*a*sa.^2)./ma + (2*sa.^4 - 2*a^2*sa.^2)./(ma.^2);
end

if length(bstab)>0
	mb	= m(bstab);
	sb	= s(bstab);
	mv(bstab)	= b - sb.^2./mb - (b*sb.^2)./(mb.^2) - (b^2*sb.^2-2*sb.^4)./(mb.^3);
	var(bstab)	= b^2 - (2*b*sb.^2)./mb + (2*sb.^4 - 2*a^2*sb.^2)./(mb.^2);
end

%keyboard
% check numerics
% keyboard
faulta	= find(mv<a);
if ~isempty(faulta)
	disp('tN:(mv<a) ');
% 	keyboard
    mv(faulta) = a;
    var(faulta) = abs(mv(faulta));
end
faultb = find(mv>b);
if ~isempty(faultb)
	disp('tN:(mv>b) ');
% 	keyboard
    mv(faultb) = b;
    var(faultb) = abs(mv(faultb));
end
end
