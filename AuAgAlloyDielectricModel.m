function [ lambda, DielectricFunction ] = AuAgAlloyDielectricModel(GMF)
% DielectricModel computes the dielectric function of a Gold-Silver alloy
% with gold molar fraction equal to GMF, according to the presented model.
% GMF takes values between 0 (pure silver) and 1 (pure gold).
% Example: [Lambda, AuAg5050] = DielectricModel(0.5);
c=2.99792458e17; % Speed of light, in nm/s
h=4.135667516e-15; % Plancks constant, in eV.s
lambda=270:1:1200; % In nm
omega=h*c./lambda; % In eV
ModelParameters=[8.9234,8.5546,9.0218;... % wp
 0.042389,0.022427,0.16713;... % gammap
 2.2715,1.7381,2.2838;... % einf
 2.6652,4.0575,3.0209;... % wg1
 2.3957,3.9260,2.7976;... % w01
 0.17880,0.017723,0.18833;... % gamma1
 73.251,51.217,22.996;... % A1
 3.5362,4.1655,3.3400;... % w02
 0.35467,0.18819,0.68309;... % gamma2
 40.007,30.770,57.540;]; % A2
wp = GMF^2*(2*ModelParameters(1,1)-4*ModelParameters(1,3)+2*ModelParameters(1,2))
+ ...
 GMF*(-ModelParameters(1,1)+4*ModelParameters(1,3)-3*ModelParameters(1,2)) + ...
 ModelParameters(1,2);

gammap = GMF^2*(2*ModelParameters(2,1)-
4*ModelParameters(2,3)+2*ModelParameters(2,2)) + ...
 GMF*(-ModelParameters(2,1)+4*ModelParameters(2,3)-3*ModelParameters(2,2)) + ...
 ModelParameters(2,2);

einf = GMF^2*(2*ModelParameters(3,1)-4*ModelParameters(3,3)+2*ModelParameters(3,2))
+ ...
 GMF*(-ModelParameters(3,1)+4*ModelParameters(3,3)-3*ModelParameters(3,2)) + ...
 ModelParameters(3,2);

wg1 = GMF^2*(2*ModelParameters(4,1)-
4*ModelParameters(4,3)+2*ModelParameters(4,2)) + ...
 GMF*(-ModelParameters(4,1)+4*ModelParameters(4,3)-3*ModelParameters(4,2)) + ...
 ModelParameters(4,2);

w01 = GMF^2*(2*ModelParameters(5,1)-
4*ModelParameters(5,3)+2*ModelParameters(5,2)) + ...
 GMF*(-ModelParameters(5,1)+4*ModelParameters(5,3)-3*ModelParameters(5,2)) + ...
 ModelParameters(5,2);

gamma1 = GMF^2*(2*ModelParameters(6,1)-
4*ModelParameters(6,3)+2*ModelParameters(6,2)) + ...
 GMF*(-ModelParameters(6,1)+4*ModelParameters(6,3)-3*ModelParameters(6,2)) + ...
 ModelParameters(6,2);

A1 = GMF^2*(2*ModelParameters(7,1)-4*ModelParameters(7,3)+2*ModelParameters(7,2))
+ ...
 GMF*(-ModelParameters(7,1)+4*ModelParameters(7,3)-3*ModelParameters(7,2)) + ...
 ModelParameters(7,2);

w02 = GMF^2*(2*ModelParameters(8,1)-
4*ModelParameters(8,3)+2*ModelParameters(8,2)) + ...
 GMF*(-ModelParameters(8,1)+4*ModelParameters(8,3)-3*ModelParameters(8,2)) + ...
 ModelParameters(8,2);

gamma2 = GMF^2*(2*ModelParameters(9,1)-
4*ModelParameters(9,3)+2*ModelParameters(9,2)) + ...
 GMF*(-ModelParameters(9,1)+4*ModelParameters(9,3)-3*ModelParameters(9,2)) + ...
 ModelParameters(9,2);

A2 = GMF^2*(2*ModelParameters(10,1)-
4*ModelParameters(10,3)+2*ModelParameters(10,2)) + ...
 GMF*(-ModelParameters(10,1)+4*ModelParameters(10,3)-3*ModelParameters(10,2)) + ...
 ModelParameters(10,2);

Drude = einf - ((wp^2)./((omega.^2) + 1i*gammap*omega));
CP1 = A1*((1./((omega+1i*gamma1).^2)) .* (-sqrt(omega+1i*gamma1-
wg1).*atan(sqrt((wg1-w01)./(omega+1i*gamma1-wg1)))) ...
 + (1./((omega+1i*gamma1).^2)) .* (-sqrt(omega+1i*gamma1+wg1).*atanh(sqrt((wg1-
w01)./(omega+1i*gamma1+wg1)))) ...
 + (1./((omega+1i*gamma1).^2)) .* (2*sqrt(wg1)*atanh(sqrt((wg1-w01)/wg1))) ...
 - sqrt(wg1-w01)*log(1-((omega+1i*gamma1)/w01).^2)./(2*(omega+1i*gamma1).^2));

CP2 = - A2 * log(1-((omega+1i*gamma2)/w02).^2)./(2*(omega+1i*gamma2).^2);

DielectricFunction = Drude + CP1 + CP2;

end