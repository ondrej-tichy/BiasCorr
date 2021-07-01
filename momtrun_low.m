function [mv,vari,nk]=tN_low(m,s,a);
% function evaluates moments of truncated normal distribution

alf	= (a-m)./(sqrt(2)*s);

% Evaluation has sense for stable only
istab	= find((alf<=3));
astab	= find(alf>3);

mv	= zeros(size(m));
vari	= mv;
nk	= mv;
%%%%%%%% STABLE
if length(istab>0)
	alf	= alf(istab);

	pom	= (1 - erf(alf))* sqrt(pi/2);
	gam	= ( - exp(-alf.^2)) ./ pom;
	del	= ( - a.*exp(-alf.^2)) ./ pom;

	mv(istab)	= m(istab) - s(istab).*gam;
	vari(istab)	= s(istab).^2 + m(istab).*mv(istab) - s(istab).*del;
	nk(istab)	= s(istab) .* pom;
end
%%%%%%%% UNSTABLE
if length(astab)>0
	ma	= m(astab);
	sa	= s(astab);
	mv(astab)	= a - sa.^2./ma;% - (a*sa.^2)./(ma.^2); % - (a^2*sa.^2-2*sa.^4)./(ma.^3);
	vari(astab)	= a^2 - (2*a*sa.^2)./ma + (2*sa.^4 - 2*a^2*sa.^2)./(ma.^2);
	nk(astab)	= -(sa.^2)./ma.*exp(-(a-ma).^2./(2*sa.^2));
end

% check numerics
fault	= find(mv<a);
if ~isempty(fault)
	disp('tN:(mv<a) ');
	keyboard
end
