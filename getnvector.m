function n = getnvector(alpha, delta)
n = [cosd(alpha); sind(alpha)*sind(delta); sind(alpha)*cosd(delta)];
end