function subsidenceGUI
% SUBSIDENCEGUI

    defaults.H = 300; defaults.drawDeg = 13; 
    defaults.NEW = 0.40; defaults.h = 4.6; defaults.e = 0.66;
    defaults.dx = 5; defaults.dy = 5;
    defaults.x0 = 100; defaults.y0 = 100;
    defaults.W = 198; defaults.L = 294;
    
    appData.leaseX = []; appData.leaseY = [];
    appData.panelX = []; appData.panelY = [];

    fig = uifigure('Name','Ground Subsidence Prediction', 'Position',[40 40 1500 850], 'Color',[0.96 0.96 0.96]);
    uilabel(fig, 'Text','GROUND SUBSIDENCE PREDICTION', 'FontSize',24, 'FontWeight','bold', 'HorizontalAlignment','center', 'Position',[400 805 900 35]);

    inputPanel = uipanel(fig, 'Title','Model Parameters', 'FontSize',14, 'Position',[20 20 330 770]);
    y = 720; dy = 30;

    makeLabel('Panel Start X (m)',y); x0field = uieditfield(inputPanel,'numeric', 'Value',defaults.x0,'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Panel Start Y (m)',y); y0field = uieditfield(inputPanel,'numeric', 'Value',defaults.y0,'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Panel Width W (m)',y); Wfield = uieditfield(inputPanel,'numeric', 'Value',defaults.W,'Limits',[0 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Panel Length L (m)',y); Lfield = uieditfield(inputPanel,'numeric', 'Value',defaults.L,'Limits',[0 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Depth H (m)',y); Hfield = uieditfield(inputPanel,'numeric', 'Value',defaults.H,'Limits',[0 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Draw Angle (deg)',y); drawField = uieditfield(inputPanel,'numeric', 'Value',defaults.drawDeg,'Limits',[0.01 89], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('NEW',y); NEWfield = uieditfield(inputPanel,'numeric', 'Value',defaults.NEW,'Limits',[0.001 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Seam thickness h (m)',y); hfield = uieditfield(inputPanel,'numeric', 'Value',defaults.h,'Limits',[0 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Extraction factor e',y); efield = uieditfield(inputPanel,'numeric', 'Value',defaults.e,'Limits',[0 1], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Grid Spacing X (m)',y); dxfield = uieditfield(inputPanel,'numeric', 'Value',defaults.dx,'Limits',[1 Inf], 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Grid Spacing Y (m)',y); dyfield = uieditfield(inputPanel,'numeric', 'Value',defaults.dy,'Limits',[1 Inf], 'Position',[175 y-4 125 25]);
    
    y=y-dy-10; uibutton(inputPanel, 'Text','1. Define Mining Panel', 'FontSize',12, 'FontWeight','bold', 'Position',[60 y 210 35], 'ButtonPushedFcn',@definePanel);
    panelStatus = uilabel(inputPanel, 'Text','Mining Panel: Using W and L fields', 'FontSize',11, 'FontColor',[0.4 0.4 0.4], 'HorizontalAlignment','center', 'Position',[20 y-25 290 25]);

    y=y-dy-20; uibutton(inputPanel, 'Text','2. Define Lease Boundary', 'FontSize',12, 'FontWeight','bold', 'Position',[60 y 210 35], 'ButtonPushedFcn',@defineBoundary);
    leaseStatus = uilabel(inputPanel, 'Text','Lease Boundary: Not Defined', 'FontSize',11, 'FontColor','red', 'HorizontalAlignment','center', 'Position',[20 y-25 290 25]);

    y=y-dy-25; makeLabel('Profile X-X at Y (m):',y); profYfield = uieditfield(inputPanel,'numeric', 'Value',150, 'Position',[175 y-4 125 25]);
    y=y-dy; makeLabel('Profile Y-Y at X (m):',y); profXfield = uieditfield(inputPanel,'numeric', 'Value',150, 'Position',[175 y-4 125 25]);

    calcButton = uibutton(inputPanel, 'Text','CALCULATE', 'FontSize',15, 'FontWeight','bold', 'Position',[60 y-50 210 40], 'Enable','off', 'ButtonPushedFcn',@calculateSubsidence);
    
    resultsArea = uitextarea(inputPanel, 'Position', [20 20 290 120], 'Editable', 'off', 'FontSize', 11);
    resultsArea.Value = {'Numerical Results Summary:','Waiting for calculation...'};

    % --- FIX: Define explicit individual axes instead of using subplot ---
    tabs = uitabgroup(fig, 'Position',[365 20 1115 770]);
    tabContour = uitab(tabs,'Title','2D Contour Map'); 
    axContour = uiaxes(tabContour,'Position',[50 65 1015 650]);
    
    tabProfX = uitab(tabs,'Title','Profile X-X (Sub, Slope, Strain)'); 
    axX_sub = uiaxes(tabProfX, 'Position', [50 490 1000 220]);
    axX_slp = uiaxes(tabProfX, 'Position', [50 270 1000 200]);
    axX_str = uiaxes(tabProfX, 'Position', [50 50 1000 200]);
    
    tabProfY = uitab(tabs,'Title','Profile Y-Y (Sub, Slope, Strain)'); 
    axY_sub = uiaxes(tabProfY, 'Position', [50 490 1000 220]);
    axY_slp = uiaxes(tabProfY, 'Position', [50 270 1000 200]);
    axY_str = uiaxes(tabProfY, 'Position', [50 50 1000 200]);
    
    tab3D = uitab(tabs,'Title','3D Subsidence'); 
    ax3D = uiaxes(tab3D,'Position',[50 65 1015 650]);

    function makeLabel(txt,yy)
        uilabel(inputPanel, 'Text',txt, 'FontSize',11, 'Position',[10 yy 160 25]);
    end

    function getPolygon(titleStr, targetX, targetY, statusObj)
        px = []; py = [];
        while true
            choice = menu([titleStr ' Input:'], 'Add Point (X,Y)', 'Finish');
            if choice == 1
                ans = inputdlg({'X (m):', 'Y (m):'}, ['Point ' num2str(length(px)+1)], [1 40], {'0','0'});
                if ~isempty(ans)
                    x_val = str2double(ans{1});
                    y_val = str2double(ans{2});
                    if isnan(x_val) || isnan(y_val)
                        uialert(fig, 'Invalid coordinate. Please enter numbers only.', 'Input Error');
                    else
                        px(end+1) = x_val; 
                        py(end+1) = y_val;
                    end
                end
            elseif choice == 2 || choice == 0
                if length(px) < 3
                    uialert(fig,'Need at least 3 valid points.','Error'); 
                else
                    break; 
                end
            end
        end
        
        if length(px) > 1
            keep = [true, (abs(diff(px)) > 1e-6) | (abs(diff(py)) > 1e-6)];
            px = px(keep); py = py(keep);
        end
        
        if px(1) ~= px(end) || py(1) ~= py(end)
            px(end+1) = px(1); py(end+1) = py(1);
        end
        
        appData.(targetX) = px; appData.(targetY) = py;
        statusObj.Text = sprintf('%s: Defined (%d pts)', titleStr, length(px));
        statusObj.FontColor = [0 0.5 0];
        
        if ~isempty(appData.leaseX), calcButton.Enable = 'on'; end
    end

    function definePanel(~,~), getPolygon('Mining Panel', 'panelX', 'panelY', panelStatus); end
    function defineBoundary(~,~), getPolygon('Lease Boundary', 'leaseX', 'leaseY', leaseStatus); end

    function calculateSubsidence(~,~)
        calcButton.Enable = 'off'; leaseStatus.Text = 'Calculating...'; drawnow;
        try
            p.W = Wfield.Value; p.L = Lfield.Value;
            p.H = Hfield.Value; p.drawDeg = drawField.Value; p.NEW = NEWfield.Value; 
            p.h = hfield.Value; p.e = efield.Value; p.dx = dxfield.Value; p.dy = dyfield.Value;
            
            if isempty(appData.panelX)
                px0 = x0field.Value; py0 = y0field.Value;
                pw = Wfield.Value; pl = Lfield.Value;
                p.panelX = [px0, px0+pw, px0+pw, px0, px0];
                p.panelY = [py0, py0, py0+pl, py0+pl, py0];
            else
                p.panelX = appData.panelX; p.panelY = appData.panelY;
            end
            p.leaseX = appData.leaseX; p.leaseY = appData.leaseY;
            
            r = subsidence_model(p);
            drawResults(r);
            updateResultsArea(r.stats);
            leaseStatus.Text = 'Status: Complete';
        catch ME
            leaseStatus.Text = 'Status: Error'; uialert(fig,ME.message,'Error');
        end
        calcButton.Enable = 'on'; drawnow;
    end

    function updateResultsArea(stats)
        resStr = {
            '--- NUMERICAL RESULTS WITHIN LEASE ---',
            sprintf('Max Subsidence: %.2f mm', stats.maxSub_pos),
            sprintf('Max Location (X,Y): (%.1f, %.1f) m', stats.maxX, stats.maxY),
            sprintf('Max Slope: %.3f mm/m', stats.maxSlope),
            sprintf('Max Tensile Strain (+): %.3f mm/m', stats.maxTensile),
            sprintf('Max Compressive Strain (-): %.3f mm/m', stats.maxCompressive)
        };
        resultsArea.Value = resStr;
    end

    function drawResults(r)
        margin = max(max(r.leaseX)-min(r.leaseX), max(r.leaseY)-min(r.leaseY)) * 0.05;
        xlims = [min(r.leaseX)-margin, max(r.leaseX)+margin];
        ylims = [min(r.leaseY)-margin, max(r.leaseY)+margin];

        Sub_mm = abs(r.Subsidence * 1000); 

        % 1. 2D Contour Map 
        cla(axContour);
        validSubs = Sub_mm(~isnan(Sub_mm));
        if isempty(validSubs) || (max(validSubs) - min(validSubs)) < 1, levels = 10; 
        else, levels = [5, 50, 100, 500, 1000, 1500, 2000, 2500]; end
        
        [C, hObj] = contour(axContour, r.X, r.Y, Sub_mm, levels, 'Color', 'k');
        clabel(C, hObj, 'FontSize', 9, 'Color', 'k');
        hold(axContour,'on');
        p_panel = plot(axContour, r.panelX, r.panelY, 'k-', 'LineWidth', 1.5);
        p_lease = plot(axContour, r.leaseX, r.leaseY, 'k-.', 'LineWidth', 1.2);
        hold(axContour,'off');
        axis(axContour, 'equal'); xlim(axContour, xlims); ylim(axContour, ylims);
        title(axContour,'(a) Predicted Subsidence Contours (mm)'); xlabel(axContour,'X, m'); ylabel(axContour,'Y, m');
        legend(axContour, [p_panel, p_lease], {'Extraction Area', 'Lease Boundary'}, 'Location', 'best');
        grid(axContour,'off');

        % 2. Multi-Profile Plot X-X
        profY = profYfield.Value;
        [~, yIdx] = min(abs(r.Y(:,1) - profY));
        profX_data = r.X(yIdx, :);
        
        cla(axX_sub); plot(axX_sub, profX_data, Sub_mm(yIdx, :), 'k.-', 'MarkerSize', 10);
        set(axX_sub, 'YDir', 'reverse'); title(axX_sub, sprintf('Profile X-X at Y = %.1f m', profY)); ylabel(axX_sub, 'Sub (mm)'); grid(axX_sub, 'on'); xlim(axX_sub, xlims);
        
        cla(axX_slp); plot(axX_slp, profX_data, r.SlopeX(yIdx, :), 'b.-', 'MarkerSize', 8);
        ylabel(axX_slp, 'Slope (mm/m)'); grid(axX_slp, 'on'); xlim(axX_slp, xlims);
        
        cla(axX_str); plot(axX_str, profX_data, r.StrainX(yIdx, :), 'r.-', 'MarkerSize', 8);
        ylabel(axX_str, 'Strain (mm/m)'); xlabel(axX_str, 'X, m'); grid(axX_str, 'on'); xlim(axX_str, xlims);

        % 3. Multi-Profile Plot Y-Y
        profX = profXfield.Value;
        [~, xIdx] = min(abs(r.X(1,:) - profX));
        profY_data = r.Y(:, xIdx);
        
        cla(axY_sub); plot(axY_sub, profY_data, Sub_mm(:, xIdx), 'k.-', 'MarkerSize', 10);
        set(axY_sub, 'YDir', 'reverse'); title(axY_sub, sprintf('Profile Y-Y at X = %.1f m', profX)); ylabel(axY_sub, 'Sub (mm)'); grid(axY_sub, 'on'); xlim(axY_sub, ylims);
        
        cla(axY_slp); plot(axY_slp, profY_data, r.SlopeY(:, xIdx), 'b.-', 'MarkerSize', 8);
        ylabel(axY_slp, 'Slope (mm/m)'); grid(axY_slp, 'on'); xlim(axY_slp, ylims);
        
        cla(axY_str); plot(axY_str, profY_data, r.StrainY(:, xIdx), 'r.-', 'MarkerSize', 8);
        ylabel(axY_str, 'Strain (mm/m)'); xlabel(axY_str, 'Y, m'); grid(axY_str, 'on'); xlim(axY_str, ylims);

        % 4. 3D Surface
        cla(ax3D);
        s = surf(ax3D, r.X, r.Y, Sub_mm); s.EdgeColor = 'none'; shading(ax3D,'interp');
        hold(ax3D,'on'); 
        plot3(ax3D, r.panelX, r.panelY, zeros(size(r.panelX)), 'k-','LineWidth',1.5); 
        hold(ax3D,'off');
        title(ax3D,'3D Subsidence Prediction (Magnitude)'); 
        xlabel(ax3D,'X (m)'); ylabel(ax3D,'Y (m)'); zlabel(ax3D,'Subsidence (mm)');
        set(ax3D, 'ZDir', 'reverse'); 
        xlim(ax3D, xlims); ylim(ax3D, ylims); grid(ax3D,'on'); view(ax3D,45,30);
    end
end