function [d_theta_ds, d_phi_ds]=angle_diff(a,b,c,theta,phi)

st=sin(theta); ct=cos(theta);
Phi_theta=(a*st)^2+(b*ct)^2;
d_theta_ds=Phi_theta^(-0.5)/sin(phi);

sp=sin(phi); cp=cos(phi);
d2=(a*ct)^2 + (b*st)^2;
Phi_phi=c^2 +d2*cp*cp;
d_phi_ds=Phi_phi^(-0.5)/sin(phi);
end