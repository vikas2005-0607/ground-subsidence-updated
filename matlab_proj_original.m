% MATLAB_PROJ: Interactive Ground Subsidence Prediction
% Separates underground mining source from surface lease boundary.

H=300; draw=(pi/180)*13; r=(H*tan(draw))/10; NEW=0.4; h=4.6; e=0.66;
W=198; L=294;

% Underground mining panel geometry
panelX = [100, 100+W, 100+W, 100];
panelY = [100, 100, 100+L, 100+L];
x_co = [panelX, panelX(1)]; y_co = [panelY, panelY(1)];

sub_factor = 0.33*(1+1.11*tanh(1.4*(W/H/NEW-1.8))); 
Svalue = [0.0297,0.0854,0.1311,0.1607,0.1708, 0.1596,0.1294,0.0859,0.0397,0.0069]/65;
Smax = sub_factor*h*e;

disp('--- Define Surface Mining Lease Boundary ---');
leaseX = []; leaseY = [];

while true
    choice = menu('Lease Area Boundary Input:', 'Add Boundary Point (X,Y)', 'Stop / Finish Input');
    if choice == 1
        answer = inputdlg({'X-coordinate (m):', 'Y-coordinate (m):'}, ...
            sprintf('Lease Point %d', length(leaseX)+1), [1 40], {'0', '0'});
        if ~isempty(answer)
            leaseX(end+1) = str2double(answer{1});
            leaseY(end+1) = str2double(answer{2});
            fprintf('Added point: (%.2f, %.2f)\n', leaseX(end), leaseY(end));
        end
    elseif choice == 2 || choice == 0
        keep = [true, (abs(diff(leaseX)) > 1e-6) | (abs(diff(leaseY)) > 1e-6)];
        leaseX = leaseX(keep); leaseY = leaseY(keep);
        if length(leaseX) < 3
            disp('Error: Polygon requires at least 3 unique points.');
        else
            break;
        end
    end
end

if leaseX(1) ~= leaseX(end) || leaseY(1) ~= leaseY(end)
    leaseX(end+1) = leaseX(1); leaseY(end+1) = leaseY(1);
end

% DYNAMIC SCALING: Calculate default grid resolution based on lease area extent
leaseWidth = max(leaseX) - min(leaseX);
leaseHeight = max(leaseY) - min(leaseY);
defaultSpacing = max(1, round(max(leaseWidth, leaseHeight) / 50));

gridAns = inputdlg({'Grid X-spacing (m):', 'Grid Y-spacing (m):'}, ...
    'Define Calculation Grid Resolution', [1 40], {num2str(defaultSpacing), num2str(defaultSpacing)});

if ~isempty(gridAns)
    grid_dx = max(1, str2double(gridAns{1}));
    grid_dy = max(1, str2double(gridAns{2}));
else
    grid_dx = defaultSpacing; grid_dy = defaultSpacing;
end

% Generate calculation grid covering strictly the bounding box of the lease area
[X, Y] = meshgrid(min(leaseX):grid_dx:max(leaseX), min(leaseY):grid_dy:max(leaseY));
inLease = inpolygon(X, Y, leaseX, leaseY);
[numRows, numCols] = size(X);
Subsidence = NaN(numRows, numCols); % NaN mask restricts calculations/plots to lease area

hWait = waitbar(0, 'Calculating Subsidence Influence...');
for m=1:numRows                         
    for n=1:numCols
        if ~inLease(m,n), continue; end
        a_ = zeros(10, 65); 
        for i=1:10                  
            theta=0:pi/32:2*pi;
            xunit=i*r*cos(theta)+X(m,n);
            yunit=i*r*sin(theta)+Y(m,n);
            for j=1:65               
                if inpolygon(xunit(j),yunit(j),panelX,panelY)
                    xx = xunit(j); yy = yunit(j);
                    db = min([xx - 100, 100+W - xx, yy - 100, 100+L - yy]);
                    d = yy - 100;
                    Wz = 0.5*tanh((5*db)/(1.5*NEW*h)-2.4)+0.5;
                    Qz = 0.9-0.1*tanh(0.5*(d-0.4*L));
                    a_(i,j) = Svalue(i)*(-Smax)*Wz*Qz*100;
                end
            end
        end
        Subsidence(m,n)=sum(sum(a_));
    end
    waitbar(m/numRows, hWait);
end
close(hWait);

% Dynamic plot limits and contour levels
margin = max(leaseWidth, leaseHeight) * 0.05;
xlims = [min(leaseX)-margin, max(leaseX)+margin];
ylims = [min(leaseY)-margin, max(leaseY)+margin];

validSubs = Subsidence(~isnan(Subsidence));
if isempty(validSubs) || (max(validSubs) - min(validSubs)) < 1e-5
    contourLevels = 15; 
else
    contourLevels = linspace(min(validSubs), max(validSubs), 15);
end

% 1. 3D Surface Plot
figure('Color', 'w');
s = surf(X, Y, Subsidence);
s.EdgeColor = 'none'; shading interp; hold on;
plot3(leaseX, leaseY, zeros(size(leaseX)), 'b-', 'LineWidth', 2);
plot3(x_co, y_co, zeros(size(x_co)), 'r--', 'LineWidth', 1.5);
xlim(xlims); ylim(ylims); 
% Note: 'axis equal' is deliberately omitted here so the Z-axis (depth) scales visibly
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Vertical Displacement / Depth (m)');
title('3D Ground Subsidence Prediction Within Lease Area');
grid on; view(45,30);

% 2. 2D Contour Plot
figure('Color', 'w');
[C, hObj] = contour(X, Y, Subsidence, contourLevels, 'Color', 'k', 'ShowText', 'on');
hold on;
p_lease = plot(leaseX, leaseY, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2);
p_panel = plot(x_co, y_co, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.5);

axis equal; % Applied here to prevent geographical distortion in 2D
xlim(xlims); ylim(ylims);

xlabel('X (m)'); ylabel('Y (m)');
title('2D Subsidence Contour Map Within Mining Lease Area');
legend([p_lease, p_panel], {'Mining Lease Boundary', 'Underground Mining Panel'}, 'Location', 'best');
grid on; set(gca, 'Color', 'w');