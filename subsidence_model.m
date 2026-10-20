function result = subsidence_model(p)
% SUBSIDENCE_MODEL
% Generalized influence-function model calculating Subsidence, Slope, and Strain fields.

    if ~isfield(p,'mode'), p.mode = "Single Seam"; end
    if ~isfield(p,'x0'), p.x0 = 100; end
    if ~isfield(p,'y0'), p.y0 = 100; end
    if ~isfield(p,'T'), p.T = 0; end
    
    if p.H <= 0 || p.NEW <= 0 || p.h <= 0 || p.e <= 0 || p.drawDeg <= 0 
        error('Mining panel dimensions and parameters must be positive.');
    end
    if ~isfield(p, 'leaseX') || ~isfield(p, 'leaseY') || length(p.leaseX) < 3
        error('A valid surface lease boundary must be provided.');
    end

    panelX = p.panelX;
    panelY = p.panelY;
    leaseX = p.leaseX; leaseY = p.leaseY;
    
    % Safety catch for NaN leakage
    valL = ~isnan(leaseX) & ~isnan(leaseY);
    leaseX = leaseX(valL); leaseY = leaseY(valL);
    valP = ~isnan(panelX) & ~isnan(panelY);
    panelX = panelX(valP); panelY = panelY(valP);

    R = p.H*tand(p.drawDeg);

    leaseWidth = max(leaseX) - min(leaseX);
    leaseHeight = max(leaseY) - min(leaseY);
    minAllowedSpacing = max([leaseWidth, leaseHeight]) / 1000;
    
    if ~isfield(p,'dx') || p.dx < minAllowedSpacing, p.dx = max(minAllowedSpacing, 1); end
    if ~isfield(p,'dy') || p.dy < minAllowedSpacing, p.dy = max(minAllowedSpacing, 1); end

    xmin = floor(min(leaseX)/p.dx)*p.dx; xmax = ceil(max(leaseX)/p.dx)*p.dx;
    ymin = floor(min(leaseY)/p.dy)*p.dy; ymax = ceil(max(leaseY)/p.dy)*p.dy;
    
    [X,Y] = meshgrid(xmin:p.dx:xmax, ymin:p.dy:ymax);
    inLease = inpolygon(X, Y, leaseX, leaseY);

    ringWeights = [0.0297, 0.0854, 0.1311, 0.1607, 0.1708, 0.1596, 0.1294, 0.0859, 0.0397, 0.0069];
    sectorWeights = ringWeights / 64;
    
    W_emp = max(panelX) - min(panelX); 
    L_emp = max(panelY) - min(panelY);
    widthRatio = W_emp/(p.H*p.NEW);

    baseFactor = 0.33 * (1 + 1.1*tanh(1.4*(widthRatio - 1.8)));
    baseFactor = max(baseFactor, 0.01); 
    Smax = baseFactor*p.h*p.e;

    Subsidence = NaN(size(X));
    Ux = NaN(size(X)); Uy = NaN(size(X));
    dr = R/10; rho = ((1:10)-0.5)*dr; theta = ((1:64)-0.5)*(2*pi/64);

    for m = 1:size(X,1)
        for n = 1:size(X,2)
            if ~inLease(m,n), continue; end 
            
            Px = X(m,n); Py = Y(m,n); Sz = 0; Sx = 0; Sy = 0;

            for i = 1:10
                xel = Px + rho(i)*cos(theta); yel = Py + rho(i)*sin(theta);
                inside = inpolygon(xel,yel,panelX,panelY);
                if ~any(inside), continue; end

                xx = xel(inside); yy = yel(inside);
                sw = sectorWeights(i);
                
                db_array = zeros(size(xx)); 
                for pt_idx = 1:length(xx)
                    pt = [xx(pt_idx), yy(pt_idx), 0];
                    num_edges = length(panelX) - 1;
                    d_edges = zeros(1, num_edges);
                    for k = 1:num_edges
                        v1 = [panelX(k), panelY(k), 0];
                        v2 = [panelX(k+1), panelY(k+1), 0];
                        d_edges(k) = point_to_line(pt, v1, v2);
                    end
                    db_array(pt_idx) = min(d_edges);
                end
                
                db = db_array;
                d = yy - min(panelY);
                
                localNEW = p.NEW*ones(size(xx));
                localSmax = Smax*ones(size(xx));

                Wz = 0.5*tanh((5*db)./(1.5*localNEW*p.h) - 2.4) + 0.5;
                Qz = 0.9 - 0.1*tanh(0.5*(d - 0.4*L_emp));
                
                dVz = -sw*localSmax.*Wz.*Qz; 
                
                Sz = Sz + sum(dVz);
                Sx = Sx + sum(dVz .* (xx - Px) / p.H); 
                Sy = Sy + sum(dVz .* (yy - Py) / p.H);
            end
            Subsidence(m,n) = Sz; 
            Ux(m,n) = Sx; 
            Uy(m,n) = Sy;
        end
    end

    [dSdy, dSdx] = gradient(Subsidence * 1000, p.dy, p.dx); 
    [dUydy, dUydx] = gradient(Uy * 1000, p.dy, p.dx);
    [dUxdy, dUxdx] = gradient(Ux * 1000, p.dy, p.dx);
    
    Exx = dUxdx; 
    Eyy = dUydy;

    validMask = ~isnan(Subsidence);
    subValid = Subsidence(validMask) * 1000; 
    slopeXValid = dSdx(validMask);
    slopeYValid = dSdy(validMask);
    strainXValid = Exx(validMask);
    strainYValid = Eyy(validMask);

    stats.maxSub = min(subValid); 
    stats.maxSub_pos = max(abs(subValid));
    
    [~, maxIdx] = min(Subsidence(:));
    [maxRow, maxCol] = ind2sub(size(Subsidence), maxIdx);
    stats.maxX = X(maxRow, maxCol);
    stats.maxY = Y(maxRow, maxCol);
    
    stats.maxSlope = max(max(abs(slopeXValid)), max(abs(slopeYValid)));
    stats.maxTensile = max(max(strainXValid), max(strainYValid));
    stats.maxCompressive = min(min(strainXValid), min(strainYValid));

    result.X = X; result.Y = Y;
    result.Subsidence = Subsidence; 
    result.SlopeX = dSdx; result.SlopeY = dSdy;
    result.StrainX = Exx; result.StrainY = Eyy;
    result.panelX = panelX; result.panelY = panelY;
    result.leaseX = leaseX; result.leaseY = leaseY;
    result.stats = stats;
end
