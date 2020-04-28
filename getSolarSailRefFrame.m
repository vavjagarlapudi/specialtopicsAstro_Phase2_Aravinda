function [r1_cap, theta_cap, eta_cap] = getSolarSailRefFrame(r1_vector)
r1_cap = r1_vector/norm(r1_vector);
theta_cap = cross([0 0 1]',r1_cap);
eta_cap = cross(r1_cap,theta_cap);
end