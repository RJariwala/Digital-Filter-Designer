classdef dsp_filter_designer < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        Panel3                        matlab.ui.container.Panel
        StopAutoButton                matlab.ui.control.Button
        StabilityNALabel_2            matlab.ui.control.Label
        Generation0Label              matlab.ui.control.Label
        TFTextArea_2                  matlab.ui.control.TextArea
        TFTextAreaLabel_2             matlab.ui.control.Label
        StopBandSlider_2              matlab.ui.control.RangeSlider
        StopBandSlider_2Label         matlab.ui.control.Label
        PassbandSlider_2              matlab.ui.control.RangeSlider
        PassbandSlider_2Label         matlab.ui.control.Label
        AutoStatusLabel               matlab.ui.control.Label
        CurrentOrderLabel             matlab.ui.control.Label
        AutoSolveButton               matlab.ui.control.Button
        PassbandRippleEditField       matlab.ui.control.NumericEditField
        PassbandRippleEditFieldLabel  matlab.ui.control.Label
        StopbandAttenEditField        matlab.ui.control.NumericEditField
        StopbandAttenEditFieldLabel   matlab.ui.control.Label
        GAORPSODropDown_2             matlab.ui.control.DropDown
        GAORPSODropDown_2Label        matlab.ui.control.Label
        FIRORIIRDropDown_2            matlab.ui.control.DropDown
        FIRORIIRDropDown_2Label       matlab.ui.control.Label
        FilterSelectDropDown_2        matlab.ui.control.DropDown
        FilterSelectDropDown_2Label   matlab.ui.control.Label
        Panel2                        matlab.ui.container.Panel
        UIAxes2                       matlab.ui.control.UIAxes
        UIAxes                        matlab.ui.control.UIAxes
        Panel                         matlab.ui.container.Panel
        MaxGenSpinner                 matlab.ui.control.Spinner
        MaxGenSpinnerLabel            matlab.ui.control.Label
        TFTextArea                    matlab.ui.control.TextArea
        TFTextAreaLabel               matlab.ui.control.Label
        StabilityNALabel              matlab.ui.control.Label
        Generation050Label            matlab.ui.control.Label
        GenerateFilterButton          matlab.ui.control.Button
        StopBandSlider                matlab.ui.control.RangeSlider
        StopBandSliderLabel           matlab.ui.control.Label
        PassbandSlider                matlab.ui.control.RangeSlider
        PassbandSliderLabel           matlab.ui.control.Label
        SIRORDERSpinner               matlab.ui.control.Spinner
        SIRORDERSpinnerLabel          matlab.ui.control.Label
        GAORPSODropDown               matlab.ui.control.DropDown
        GAORPSODropDownLabel          matlab.ui.control.Label
        FIRORIIRDropDown              matlab.ui.control.DropDown
        FIRORIIRDropDownLabel         matlab.ui.control.Label
        FilterSelectDropDown          matlab.ui.control.DropDown
        FilterSelectDropDownLabel     matlab.ui.control.Label
    end

    
    properties (Access = private)
        K_points = 512          % Number of frequency evaluation points (0 to pi)
        OmegaGrid               % Normalized frequency vector [0, pi]
        IsOptimizing = false    % Guard flag to avoid overlapping runs
        StopRequested = false   % <--- ADD THIS LINE: Interrupt flag for Stop button
        UnitCircleX             % Unit circle x-data
        UnitCircleY             % Unit circle y-data
    end
    
    methods (Access = private)
        % Live visual update routine for Panel 3
        function updateVisuals_2(app, b, a, isIIR, currentOrder, gen, maxGen, targetAs)
            % 1. Update Generation Label on Panel 3
            app.Generation0Label.Text = sprintf('Order N = %d | Gen: %d', currentOrder, gen);

            % 2. Update Stability Label on Panel 3
            if isIIR
                poles = roots(a);
                if isempty(poles)
                    maxR = 0;
                else
                    maxR = max(abs(poles));
                end
                if maxR < 1.0
                    app.StabilityNALabel_2.Text = sprintf('Stability: STABLE (|p|max = %.2f)', maxR);
                    app.StabilityNALabel_2.FontColor = [0.1 0.85 0.1];
                else
                    app.StabilityNALabel_2.Text = sprintf('Stability: UNSTABLE (|p|max = %.2f)', maxR);
                    app.StabilityNALabel_2.FontColor = [1.0 0.2 0.2];
                end
            else
                poles = zeros(length(b)-1, 1);
                app.StabilityNALabel_2.Text = 'Stability: FIR ALWAYS STABLE';
                app.StabilityNALabel_2.FontColor = [0.1 0.85 0.1];
            end

            % 3. Plot Frequency Response (dB)
            w_norm = app.OmegaGrid / pi;
            H = abs(app.evaluateFreqz(b, a, app.OmegaGrid));
            H_dB = 20 * log10(max(H, 1e-4));

            plot(app.UIAxes, w_norm, H_dB, 'c-', 'LineWidth', 1.8);
            hold(app.UIAxes, 'on');
            yline(app.UIAxes, 0, 'g:', 'LineWidth', 1.1);
            yline(app.UIAxes, -targetAs, 'r--', 'LineWidth', 1.1);
            hold(app.UIAxes, 'off');

            xlim(app.UIAxes, [0 1]);
            ylim(app.UIAxes, [-80 5]);
            grid(app.UIAxes, 'on');

            % 4. Plot Pole-Zero Map (Z-Plane)
            z_zeros = roots(b);
            plot(app.UIAxes2, app.UnitCircleX, app.UnitCircleY, 'w--', 'LineWidth', 1.1);
            hold(app.UIAxes2, 'on');
            xline(app.UIAxes2, 0, 'Color', [0.4 0.4 0.4]);
            yline(app.UIAxes2, 0, 'Color', [0.4 0.4 0.4]);

            if ~isempty(z_zeros)
                plot(app.UIAxes2, real(z_zeros), imag(z_zeros), 'co', 'MarkerSize', 7, 'LineWidth', 1.5);
            end
            if ~isempty(poles)
                plot(app.UIAxes2, real(poles), imag(poles), 'rx', 'MarkerSize', 8, 'LineWidth', 2.0);
            end
            hold(app.UIAxes2, 'off');

            axis(app.UIAxes2, 'equal');
            xlim(app.UIAxes2, [-1.5 1.5]);
            ylim(app.UIAxes2, [-1.5 1.5]);
            grid(app.UIAxes2, 'on');

            drawnow;
        end
        % Generate binary evaluation masks using Panel 3 sliders
        function [passMask, stopMask] = getFrequencyMasks_2(app)
            w_norm = app.OmegaGrid / pi;
            fType = upper(strtrim(app.FilterSelectDropDown_2.Value));
            pb = app.PassbandSlider_2.Value;
            sb = app.StopBandSlider_2.Value;

            passMask = false(size(w_norm));
            stopMask = false(size(w_norm));

            switch fType
                case 'LOW PASS'
                    passMask = (w_norm <= pb(2));
                    stopMask = (w_norm >= sb(1));
                case 'HIGH PASS'
                    stopMask = (w_norm <= sb(2));
                    passMask = (w_norm >= pb(1));
                case 'BAND PASS'
                    passMask = (w_norm >= pb(1)) & (w_norm <= pb(2));
                    stopMask = (w_norm <= sb(1)) | (w_norm >= sb(2));
                case 'BAND STOP'
                    passMask = (w_norm <= pb(1)) | (w_norm >= pb(2));
                    stopMask = (w_norm >= sb(1)) & (w_norm <= sb(2));
            end
        end

        % Estimate baseline order using Panel 3 slider transition width
        function Nest = estimateInitialOrder_2(app, As_dB, Rp_dB, isIIR)
            pb = app.PassbandSlider_2.Value;
            sb = app.StopBandSlider_2.Value;
            fType = upper(strtrim(app.FilterSelectDropDown_2.Value));

            switch fType
                case 'LOW PASS'
                    df = abs(sb(1) - pb(2)) / 2;
                case 'HIGH PASS'
                    df = abs(pb(1) - sb(2)) / 2;
                case {'BAND PASS', 'BAND STOP'}
                    df = min(abs(pb(1) - sb(1)), abs(sb(2) - pb(2))) / 2;
            end
            df = max(0.02, df);

            if ~isIIR
                % Bellanger equation for FIR
                delta_p = (10^(Rp_dB / 20) - 1) / (10^(Rp_dB / 20) + 1);
                delta_s = 10^(-As_dB / 20);
                Nest = ceil((-20 * log10(sqrt(delta_p * delta_s)) - 13) / (14.6 * df));
                Nest = max(2, min(Nest, 24));
                if mod(Nest, 2) ~= 0, Nest = Nest + 1; end % Enforce even order for FIR
            else
                % Analytical approximation for IIR
                Nest = ceil(As_dB / (14.0 * log10(1 + df * 10)));
                Nest = max(2, min(Nest, 10));
                if mod(Nest, 2) ~= 0, Nest = Nest + 1; end
            end
        end
        % Enforce safe slider boundaries for Panel 3
        function updateSliderEdges_2(app)
            fType = upper(strtrim(app.FilterSelectDropDown_2.Value));
            pb = app.PassbandSlider_2.Value;
            sb = app.StopBandSlider_2.Value;

            switch fType
                case 'LOW PASS'
                    wp = pb(2);
                    if wp <= 0.02 || wp >= 0.95, wp = 0.35; end
                    ws = max(wp + 0.05, sb(1));
                    ws = min(max(ws, wp + 0.02), 0.98);
                    app.PassbandSlider_2.Value = [0, wp];
                    app.StopBandSlider_2.Value = [ws, 1];

                case 'HIGH PASS'
                    ws = sb(2);
                    if ws <= 0.02 || ws >= 0.95, ws = 0.35; end
                    wp = max(ws + 0.05, pb(1));
                    wp = min(max(wp, ws + 0.02), 0.98);
                    app.StopBandSlider_2.Value = [0, ws];
                    app.PassbandSlider_2.Value = [wp, 1];

                case 'BAND PASS'
                    wp1 = pb(1); wp2 = pb(2);
                    ws1 = sb(1); ws2 = sb(2);
                    if wp1 <= 0.01 || wp2 >= 0.99 || ws1 <= 0.01 || ws2 >= 0.99 || wp1 >= wp2 || ws1 >= ws2
                        wp1 = 0.35; wp2 = 0.65;
                        ws1 = 0.20; ws2 = 0.80;
                    else
                        wp1 = max(0.08, min(wp1, 0.80));
                        wp2 = max(wp1 + 0.05, min(wp2, 0.90));
                        ws1 = max(0.02, min(ws1, wp1 - 0.03));
                        ws2 = min(0.98, max(ws2, wp2 + 0.03));
                    end
                    app.PassbandSlider_2.Value = [wp1, wp2];
                    app.StopBandSlider_2.Value = [ws1, ws2];

                case 'BAND STOP'
                    wp1 = pb(1); wp2 = pb(2);
                    ws1 = sb(1); ws2 = sb(2);
                    if wp1 <= 0.01 || wp2 >= 0.99 || ws1 <= 0.01 || ws2 >= 0.99 || wp1 >= wp2 || ws1 >= ws2
                        ws1 = 0.35; ws2 = 0.65;
                        wp1 = 0.20; wp2 = 0.80;
                    else
                        ws1 = max(0.08, min(ws1, 0.80));
                        ws2 = max(ws1 + 0.05, min(ws2, 0.90));
                        wp1 = max(0.02, min(wp1, ws1 - 0.03));
                        wp2 = min(0.98, max(ws2, wp2 + 0.03));
                    end
                    app.PassbandSlider_2.Value = [wp1, wp2];
                    app.StopBandSlider_2.Value = [ws1, ws2];
            end
        end
        % Estimates minimum baseline filter order from user specs and band edges
        function Nest = estimateInitialOrder(app, As_dB, Rp_dB, isIIR, fType)
            pb = app.PassbandSlider.Value;
            sb = app.StopBandSlider.Value;

            switch fType
                case 'LOW PASS'
                    df = abs(sb(1) - pb(2)) / 2;
                case 'HIGH PASS'
                    df = abs(pb(1) - sb(2)) / 2;
                case {'BAND PASS', 'BAND STOP'}
                    df = min(abs(pb(1) - sb(1)), abs(sb(2) - pb(2))) / 2;
            end
            df = max(0.02, df);

            if ~isIIR
                % Bellanger approximation for FIR
                delta_p = (10^(Rp_dB / 20) - 1) / (10^(Rp_dB / 20) + 1);
                delta_s = 10^(-As_dB / 20);
                Nest = ceil((-20 * log10(sqrt(delta_p * delta_s)) - 13) / (14.6 * df));
                Nest = max(2, min(Nest, 24));
                if mod(Nest, 2) ~= 0, Nest = Nest + 1; end % Enforce even order for FIR
            else
                % Typical IIR prototype approximation
                Nest = ceil(As_dB / (14.0 * log10(1 + df * 10)));
                Nest = max(2, min(Nest, 10));
                if mod(Nest, 2) ~= 0, Nest = Nest + 1; end
            end
        end
        % Gets passband and stopband masks for a specified filter type string
        function [passMask, stopMask] = getFrequencyMasksForType(app, fType)
            w_norm = app.OmegaGrid / pi;
            pb = app.PassbandSlider.Value;
            sb = app.StopBandSlider.Value;
            passMask = false(size(w_norm));
            stopMask = false(size(w_norm));

            switch fType
                case 'LOW PASS'
                    passMask = (w_norm <= pb(2));
                    stopMask = (w_norm >= sb(1));
                case 'HIGH PASS'
                    stopMask = (w_norm <= sb(2));
                    passMask = (w_norm >= pb(1));
                case 'BAND PASS'
                    passMask = (w_norm >= pb(1)) & (w_norm <= pb(2));
                    stopMask = (w_norm <= sb(1)) | (w_norm >= sb(2));
                case 'BAND STOP'
                    passMask = (w_norm <= pb(1)) | (w_norm >= pb(2));
                    stopMask = (w_norm >= sb(1)) & (w_norm <= sb(2));
            end
        end
        % Formats coefficient vectors b and a into a clean H(z) string
        function str = formatTransferFunction(~, b, a)
            numStr = formatPoly(b);
            if length(a) == 1 && a(1) == 1
                % FIR Filter
                str = sprintf('H(z) = %s', numStr);
            else
                % IIR Filter
                denStr = formatPoly(a);
                separator = repmat('-', 1, max(length(numStr), length(denStr)) + 4);
                str = sprintf('       %s\nH(z) = %s\n       %s', numStr, separator, denStr);
            end

            function pStr = formatPoly(coeffs)
                terms = {};
                for k = 1:length(coeffs)
                    c = coeffs(k);
                    if abs(c) < 1e-4, continue; end
                    order = k - 1;
                    if order == 0
                        term = sprintf('%.4f', c);
                    elseif order == 1
                        if c >= 0
                            term = sprintf('+ %.4f z^-1', c);
                        else
                            term = sprintf('- %.4f z^-1', abs(c));
                        end
                    else
                        if c >= 0
                            term = sprintf('+ %.4f z^-%d', c, order);
                        else
                            term = sprintf('- %.4f z^-%d', abs(c), order);
                        end
                    end
                    terms{end+1} = term; %#ok<AGROW>
                end
                if isempty(terms)
                    pStr = '0';
                else
                    pStr = strtrim(strjoin(terms, ' '));
                end
            end
        end
        % Custom standalone replacement for freqz (requires no Signal Processing Toolbox)
        function H = evaluateFreqz(~, b, a, w)
            % Ensure column vectors for coefficients
            b = b(:);
            a = a(:);
            w = w(:)'; % Row vector for frequencies

            % Exponent grids: e^(-j * k * w)
            kb = (0 : length(b)-1)';
            ka = (0 : length(a)-1)';

            % Compute DTFT numerator and denominator sums via matrix multiplication
            num = b.' * exp(-1j * (kb * w));
            den = a.' * exp(-1j * (ka * w));

            H = num ./ den;
        end
        
        % Enforce RangeSlider thumb locking based on selected filter type
        % Enforce safe RangeSlider edges without violating [v1 <= v2] or limits [0, 1]
        function updateSliderEdges(app)
            fType = upper(strtrim(app.FilterSelectDropDown.Value));
            pb = app.PassbandSlider.Value;
            sb = app.StopBandSlider.Value;

            switch fType
                case 'LOW PASS'
                    % Passband: [0, wp], Stopband: [ws, 1]
                    wp = pb(2);
                    if wp <= 0.02 || wp >= 0.95
                        wp = 0.35; % Safe default
                    end
                    ws = max(wp + 0.05, sb(1));
                    ws = min(max(ws, wp + 0.02), 0.98);

                    app.PassbandSlider.Value = [0, wp];
                    app.StopBandSlider.Value = [ws, 1];

                case 'HIGH PASS'
                    % Stopband: [0, ws], Passband: [wp, 1]
                    ws = sb(2);
                    if ws <= 0.02 || ws >= 0.95
                        ws = 0.35; % Safe default
                    end
                    wp = max(ws + 0.05, pb(1));
                    wp = min(max(wp, ws + 0.02), 0.98);

                    app.StopBandSlider.Value = [0, ws];
                    app.PassbandSlider.Value = [wp, 1];

                case 'BAND PASS'
                    % Passband: [wp1, wp2], Stopband edges: [ws1, ws2]
                    % Safe default ordering: 0 < ws1 < wp1 < wp2 < ws2 < 1
                    wp1 = pb(1); wp2 = pb(2);
                    ws1 = sb(1); ws2 = sb(2);

                    % If coming from LP/HP, initialize to a clean preset
                    if wp1 <= 0.01 || wp2 >= 0.99 || ws1 <= 0.01 || ws2 >= 0.99 || wp1 >= wp2 || ws1 >= ws2
                        wp1 = 0.35; wp2 = 0.65;
                        ws1 = 0.20; ws2 = 0.80;
                    else
                        % Ensure strict ordering: ws1 < wp1 < wp2 < ws2
                        wp1 = max(0.08, min(wp1, 0.80));
                        wp2 = max(wp1 + 0.05, min(wp2, 0.90));
                        ws1 = max(0.02, min(ws1, wp1 - 0.03));
                        ws2 = min(0.98, max(ws2, wp2 + 0.03));
                    end

                    app.PassbandSlider.Value = [wp1, wp2];
                    app.StopBandSlider.Value = [ws1, ws2];

                case 'BAND STOP'
                    % Stopband: [ws1, ws2], Passband edges: [wp1, wp2]
                    % Safe default ordering: 0 < wp1 < ws1 < ws2 < wp2 < 1
                    wp1 = pb(1); wp2 = pb(2);
                    ws1 = sb(1); ws2 = sb(2);

                    % If coming from another type, initialize to a clean preset
                    if wp1 <= 0.01 || wp2 >= 0.99 || ws1 <= 0.01 || ws2 >= 0.99 || wp1 >= wp2 || ws1 >= ws2
                        ws1 = 0.35; ws2 = 0.65;
                        wp1 = 0.20; wp2 = 0.80;
                    else
                        % Ensure strict ordering: wp1 < ws1 < ws2 < wp2
                        ws1 = max(0.08, min(ws1, 0.80));
                        ws2 = max(ws1 + 0.05, min(ws2, 0.90));
                        wp1 = max(0.02, min(wp1, ws1 - 0.03));
                        wp2 = min(0.98, max(wp2, ws2 + 0.03));
                    end

                    app.PassbandSlider.Value = [wp1, wp2];
                    app.StopBandSlider.Value = [ws1, ws2];
            end
        end

        % Build binary evaluation masks for Passband and Stopband
        function [passMask, stopMask] = getFrequencyMasks(app)
            w_norm = app.OmegaGrid / pi;
            fType = upper(strtrim(app.FilterSelectDropDown.Value));
            pb = app.PassbandSlider.Value;
            sb = app.StopBandSlider.Value;

            passMask = false(size(w_norm));
            stopMask = false(size(w_norm));

            switch fType
                case 'LOW PASS'
                    passMask = (w_norm <= pb(2));
                    stopMask = (w_norm >= sb(1));
                case 'HIGH PASS'
                    stopMask = (w_norm <= sb(2));
                    passMask = (w_norm >= pb(1));
                case 'BAND PASS'
                    passMask = (w_norm >= pb(1)) & (w_norm <= pb(2));
                    stopMask = (w_norm <= sb(1)) | (w_norm >= sb(2));
                case 'BAND STOP'
                    passMask = (w_norm <= pb(1)) | (w_norm >= pb(2));
                    stopMask = (w_norm >= sb(1)) & (w_norm <= sb(2));
            end
        end

        % Cost / Fitness Function
        % Fitness Function with dB-scaled stopband penalty
        % Cost / Fitness Function
        % Standard DSP Fitness Function: Mean Squared Error against Ideal Mask
        % Standard DSP Fitness Function: Mean Squared Error against Ideal Mask
        function cost = computeFitness(app, x, N, isIIR, passMask, stopMask)
            if ~isIIR
                % FIR: Coefficients are x directly
                b = x;
                a = 1;

                % Do not normalize by sum(b).
                % DC normalization is incorrect for high-pass
                % and band-pass filters.
            else
                % IIR: Numerator and Denominator
                b = x(1 : N+1);
                a = [1, x(N+2 : end)];

                % Stability guarantee
                p_roots = roots(a);
                r = abs(p_roots);

                if any(r >= 0.99)
                    cost = 1e5 + 1e3 * sum(max(0, r - 0.95).^2);
                    return;
                end
            end

            % Compute frequency response magnitude
            H = abs(app.evaluateFreqz(b, a, app.OmegaGrid));

            % 1. Passband error: Drive |H| toward 1.0
            pbH = H(passMask);
            if isempty(pbH)
                passCost = 0;
            else
                passCost = mean((pbH - 1.0).^2);
            end

            % 2. Stopband error: Drive |H| toward 0.0
            sbH = H(stopMask);
            if isempty(sbH)
                stopCost = 0;
            else
                stopCost = mean(sbH.^2) * 2.0; % Slightly higher weight for stopband rejection
            end

            cost = passCost + stopCost;
        end

        % Live update routine called each generation
        % Live Visual Update per Generation
        function updateVisuals(app, b, a, isIIR, gen, maxGen)
            % 1. Update Generation Label
            app.Generation050Label.Text = sprintf('Generation: %d / %d', gen, maxGen);

            % 2. Stability Status
            if isIIR
                poles = roots(a);
                if isempty(poles)
                    maxR = 0;
                else
                    maxR = max(abs(poles));
                end

                if maxR < 1.0
                    app.StabilityNALabel.Text = sprintf('Stability: STABLE (|p|max = %.2f)', maxR);
                    app.StabilityNALabel.FontColor = [0.1 0.85 0.1]; % Green
                else
                    app.StabilityNALabel.Text = sprintf('Stability: UNSTABLE (|p|max = %.2f)', maxR);
                    app.StabilityNALabel.FontColor = [1.0 0.2 0.2];  % Red
                end
            else
                poles = zeros(length(b)-1, 1);
                app.StabilityNALabel.Text = 'Stability: FIR ALWAYS STABLE';
                app.StabilityNALabel.FontColor = [0.1 0.85 0.1];
            end

            % 3. Plot Frequency Response (dB)
            w_norm = app.OmegaGrid / pi;
            H = abs(app.evaluateFreqz(b, a, app.OmegaGrid));
            H_dB = 20 * log10(max(H, 1e-4)); % Floor at -80 dB

            plot(app.UIAxes, w_norm, H_dB, 'c-', 'LineWidth', 1.8);
            hold(app.UIAxes, 'on');
            yline(app.UIAxes, 0, 'g:', 'LineWidth', 1.1); % Ideal 0 dB passband reference
            yline(app.UIAxes, -app.StopbandAttenEditField.Value, 'r--', 'LineWidth', 1.1);
            hold(app.UIAxes, 'off');

            xlim(app.UIAxes, [0 1]);
            ylim(app.UIAxes, [-80 5]);
            grid(app.UIAxes, 'on');

            % 4. Plot Pole-Zero Map (Z-Plane)
            z_zeros = roots(b);

            plot(app.UIAxes2, app.UnitCircleX, app.UnitCircleY, 'w--', 'LineWidth', 1.1);
            hold(app.UIAxes2, 'on');
            xline(app.UIAxes2, 0, 'Color', [0.4 0.4 0.4]);
            yline(app.UIAxes2, 0, 'Color', [0.4 0.4 0.4]);

            % Plot Zeros (circles)
            if ~isempty(z_zeros)
                plot(app.UIAxes2, real(z_zeros), imag(z_zeros), 'co', 'MarkerSize', 7, 'LineWidth', 1.5);
            end

            % Plot Poles (crosses)
            if ~isempty(poles)
                plot(app.UIAxes2, real(poles), imag(poles), 'rx', 'MarkerSize', 8, 'LineWidth', 2.0);
            end
            hold(app.UIAxes2, 'off');

            axis(app.UIAxes2, 'equal');
            xlim(app.UIAxes2, [-1.5 1.5]);
            ylim(app.UIAxes2, [-1.5 1.5]);
            grid(app.UIAxes2, 'on');

            drawnow;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.OmegaGrid = linspace(0, pi, app.K_points);
            theta = linspace(0, 2*pi, 300);
            app.UnitCircleX = cos(theta);
            app.UnitCircleY = sin(theta);
            
            updateSliderEdges(app);
            updateSliderEdges_2(app);
            
        end

        % Value changed function: FilterSelectDropDown
        function FilterSelectDropDownValueChanged(app, event)
            value = app.FilterSelectDropDown.Value;
            updateSliderEdges(app);
        end

        % Value changed function: PassbandSlider
        function PassbandSliderValueChanged(app, event)
            value = app.PassbandSlider.Value;
            updateSliderEdges(app);
        end

        % Value changed function: StopBandSlider
        function StopBandSliderValueChanged(app, event)
            value = app.StopBandSlider.Value;
            updateSliderEdges(app);
        end

        % Button pushed function: GenerateFilterButton
        function GenerateFilterButtonPushed(app, event)
            
            if app.IsOptimizing
                return;
            end
            app.IsOptimizing = true;
            app.GenerateFilterButton.Enable = 'off';

            % 1. Read Inputs
            N = double(app.SIRORDERSpinner.Value);
            isIIR = strcmp(upper(strtrim(app.FIRORIIRDropDown.Value)), 'IIR');
            isPSO = contains(upper(app.GAORPSODropDown.Value), 'PSO');

            [passMask, stopMask] = getFrequencyMasks(app);

            % 2. Setup Variables
            if ~isIIR
                nvars = N + 1;
                lb = -2.0 * ones(1, nvars);
                ub =  2.0 * ones(1, nvars);
            else
                nvars = 2*N + 1;
                lb = [-2.0*ones(1, N+1), -1.2*ones(1, N)];
                ub = [ 2.0*ones(1, N+1),  1.2*ones(1, N)];
            end

            maxGen = double(app.MaxGenSpinner.Value);
            popSize = 40;

            % 3. Run Selected Algorithm
            if isPSO
                % Particle Swarm Optimization
                w = 0.72; c1 = 1.49; c2 = 1.49;
                particles = repmat(lb, popSize, 1) + rand(popSize, nvars) .* repmat(ub - lb, popSize, 1);
                velocities = zeros(popSize, nvars);
                pBest = particles;
                pBestScores = inf(popSize, 1);

                for i = 1:popSize
                    pBestScores(i) = computeFitness(app, particles(i,:), N, isIIR, passMask, stopMask);
                end
                [gBestScore, bestIdx] = min(pBestScores);
                gBest = pBest(bestIdx, :);

                for gen = 1:maxGen
                    for i = 1:popSize
                        r1 = rand(1, nvars); r2 = rand(1, nvars);
                        velocities(i,:) = w * velocities(i,:) ...
                            + c1 * r1 .* (pBest(i,:) - particles(i,:)) ...
                            + c2 * r2 .* (gBest - particles(i,:));
                        particles(i,:) = particles(i,:) + velocities(i,:);
                        particles(i,:) = max(min(particles(i,:), ub), lb);

                        score = computeFitness(app, particles(i,:), N, isIIR, passMask, stopMask);
                        if score < pBestScores(i)
                            pBest(i,:) = particles(i,:);
                            pBestScores(i) = score;
                            if score < gBestScore
                                gBest = particles(i,:);
                                gBestScore = score;
                            end
                        end
                    end

                    if ~isIIR
                        b_curr = gBest;
                        a_curr = 1;
                    else
                        b_curr = gBest(1:N+1);
                        a_curr = [1, gBest(N+2:end)];
                    end

                    updateVisuals(app, b_curr, a_curr, isIIR, gen, maxGen);
                end
            else
                % Genetic Algorithm (GA)
                pop = repmat(lb, popSize, 1) + rand(popSize, nvars) .* repmat(ub - lb, popSize, 1);
                scores = inf(popSize, 1);

                for i = 1:popSize
                    scores(i) = computeFitness(app, pop(i,:), N, isIIR, passMask, stopMask);
                end

                for gen = 1:maxGen
                    [scores, sortIdx] = sort(scores, 'ascend');
                    pop = pop(sortIdx, :);

                    newPop = zeros(size(pop));
                    newPop(1:2, :) = pop(1:2, :); % Elitism (keep top 2)

                    for i = 3:2:popSize
                        idx1 = min(randperm(popSize, 3));
                        idx2 = min(randperm(popSize, 3));
                        p1 = pop(idx1, :); p2 = pop(idx2, :);

                        alpha = rand(1, nvars);
                        c1 = alpha .* p1 + (1 - alpha) .* p2;
                        c2 = (1 - alpha) .* p1 + alpha .* p2;

                        if rand < 0.35, c1 = c1 + 0.12 * randn(1, nvars); end
                        if rand < 0.35, c2 = c2 + 0.12 * randn(1, nvars); end

                        newPop(i, :) = max(min(c1, ub), lb);
                        if i+1 <= popSize
                            newPop(i+1, :) = max(min(c2, ub), lb);
                        end
                    end

                    pop = newPop;
                    for i = 1:popSize
                        scores(i) = computeFitness(app, pop(i,:), N, isIIR, passMask, stopMask);
                    end

                    [~, bestIdx] = min(scores);
                    bestInd = pop(bestIdx, :);

                    if ~isIIR
                        b_curr = bestInd;
                        a_curr = 1;
                    else
                        b_curr = bestInd(1:N+1);
                        a_curr = [1, bestInd(N+2:end)];
                    end

                    updateVisuals(app, b_curr, a_curr, isIIR, gen, maxGen);
                end
            end

            % Display final Transfer Function
            tfString = app.formatTransferFunction(b_curr, a_curr);
            app.TFTextArea.Value = tfString;

            app.GenerateFilterButton.Enable = 'on';
            app.IsOptimizing = false;
        end

        % Button pushed function: AutoSolveButton
        function AutoSolveButtonPushed(app, event)

            % Button pushed function: AutoSolveButton (Panel 3 Generate Filter)
      
            
            if app.IsOptimizing
                return;
            end
            app.IsOptimizing = true;
            app.StopRequested = false;

            % UI State: Disable run buttons, enable Stop button
            app.GenerateFilterButton.Enable = 'off';
            app.AutoSolveButton.Enable = 'off';
            app.StopAutoButton.Enable = 'on';

            % 1. Read Inputs from Panel 3
            isIIR = strcmp(upper(strtrim(app.FIRORIIRDropDown_2.Value)), 'IIR');
            isPSO = contains(upper(app.GAORPSODropDown_2.Value), 'PSO');
            targetAs = max(5.0, app.StopbandAttenEditField.Value);
            targetRp = max(0.05, app.PassbandRippleEditField.Value);

            [passMask, stopMask] = app.getFrequencyMasks_2();

            % 2. Deterministically Compute the Exact Required Order
            N = app.estimateInitialOrder_2(targetAs, targetRp, isIIR);
            app.CurrentOrderLabel.Text = sprintf('Determined Order: N = %d', N);
            app.AutoStatusLabel.Text = 'Status: Optimizing toward specs...';
            app.AutoStatusLabel.FontColor = [0.9 0.8 0.1]; % Yellow

            % 3. Set Optimizer Bounds for Determined Order N
            if ~isIIR
                nvars = N + 1;
                lb = -2.0 * ones(1, nvars);
                ub =  2.0 * ones(1, nvars);
            else
                nvars = 2 * N + 1;
                lb = [-2.0 * ones(1, N + 1), -1.2 * ones(1, N)];
                ub = [ 2.0 * ones(1, N + 1),  1.2 * ones(1, N)];
            end

            maxGenLimit = 2000; % High upper ceiling so user can stop anytime
            popSize = 50;
            specMet = false;
            b_curr = [];
            a_curr = [];

            if isPSO
                % --- Continuous PSO Engine ---
                w = 0.72; c1 = 1.49; c2 = 1.49;
                particles = repmat(lb, popSize, 1) + rand(popSize, nvars) .* repmat(ub - lb, popSize, 1);
                velocities = zeros(popSize, nvars);
                pBest = particles;
                pBestScores = inf(popSize, 1);

                for i = 1:popSize
                    pBestScores(i) = computeFitness(app, particles(i,:), N, isIIR, passMask, stopMask);
                end
                [gBestScore, bestIdx] = min(pBestScores);
                gBest = pBest(bestIdx, :);

                for gen = 1:maxGenLimit
                    % Check user stop interrupt
                    if app.StopRequested
                        break;
                    end

                    for i = 1:popSize
                        r1 = rand(1, nvars); r2 = rand(1, nvars);
                        velocities(i,:) = w * velocities(i,:) ...
                            + c1 * r1 .* (pBest(i,:) - particles(i,:)) ...
                            + c2 * r2 .* (gBest - particles(i,:));
                        particles(i,:) = particles(i,:) + velocities(i,:);
                        particles(i,:) = max(min(particles(i,:), ub), lb);

                        score = computeFitness(app, particles(i,:), N, isIIR, passMask, stopMask);
                        if score < pBestScores(i)
                            pBest(i,:) = particles(i,:);
                            pBestScores(i) = score;
                            if score < gBestScore
                                gBest = particles(i,:);
                                gBestScore = score;
                            end
                        end
                    end

                    if ~isIIR
                        b_curr = gBest;
                        H_test = abs(app.evaluateFreqz(b_curr, 1, app.OmegaGrid));
                        maxPB = max(H_test(passMask));
                        if maxPB > 1e-4, b_curr = b_curr / maxPB; end
                        a_curr = 1;
                    else
                        b_curr = gBest(1:N+1);
                        a_curr = [1, gBest(N+2:end)];
                    end

                    % Live UI Update
                    updateVisuals_2(app, b_curr, a_curr, isIIR, N, gen, maxGenLimit, targetAs);

                    % Spec verification check in dB
                    H = abs(app.evaluateFreqz(b_curr, a_curr, app.OmegaGrid));
                    H_dB = 20 * log10(max(H, 1e-4));
                    actualStopMax = max(H_dB(stopMask));
                    actualPassRipple = max(abs(H_dB(passMask)));

                    if (actualStopMax <= -targetAs) && (actualPassRipple <= targetRp)
                        specMet = true;
                        app.AutoStatusLabel.Text = sprintf('Specs Met at Gen %d!', gen);
                        app.AutoStatusLabel.FontColor = [0.1 0.9 0.1]; % Green
                        break;
                    end
                end

            else
                % --- Continuous GA Engine ---
                pop = repmat(lb, popSize, 1) + rand(popSize, nvars) .* repmat(ub - lb, popSize, 1);
                scores = inf(popSize, 1);
                for i = 1:popSize
                    scores(i) = computeFitness(app, pop(i,:), N, isIIR, passMask, stopMask);
                end

                for gen = 1:maxGenLimit
                    if app.StopRequested
                        break;
                    end

                    [scores, sortIdx] = sort(scores, 'ascend');
                    pop = pop(sortIdx, :);
                    newPop = zeros(size(pop));
                    newPop(1:2, :) = pop(1:2, :);

                    for i = 3:2:popSize
                        idx1 = min(randperm(popSize, 3));
                        idx2 = min(randperm(popSize, 3));
                        p1 = pop(idx1, :); p2 = pop(idx2, :);
                        alpha = rand(1, nvars);
                        c1 = alpha .* p1 + (1 - alpha) .* p2;
                        c2 = (1 - alpha) .* p1 + alpha .* p2;
                        if rand < 0.35, c1 = c1 + 0.10 * randn(1, nvars); end
                        if rand < 0.35, c2 = c2 + 0.10 * randn(1, nvars); end
                        newPop(i, :) = max(min(c1, ub), lb);
                        if i+1 <= popSize
                            newPop(i+1, :) = max(min(c2, ub), lb);
                        end
                    end
                    pop = newPop;

                    for i = 1:popSize
                        scores(i) = computeFitness(app, pop(i,:), N, isIIR, passMask, stopMask);
                    end

                    [~, bestIdx] = min(scores);
                    bestInd = pop(bestIdx, :);

                    if ~isIIR
                        b_curr = bestInd;
                        H_test = abs(app.evaluateFreqz(b_curr, 1, app.OmegaGrid));
                        maxPB = max(H_test(passMask));
                        if maxPB > 1e-4, b_curr = b_curr / maxPB; end
                        a_curr = 1;
                    else
                        b_curr = bestInd(1:N+1);
                        a_curr = [1, bestInd(N+2:end)];
                    end

                    % Live UI Update
                    updateVisuals_2(app, b_curr, a_curr, isIIR, N, gen, maxGenLimit, targetAs);

                    % Spec verification check in dB
                    H = abs(app.evaluateFreqz(b_curr, a_curr, app.OmegaGrid));
                    H_dB = 20 * log10(max(H, 1e-4));
                    actualStopMax = max(H_dB(stopMask));
                    actualPassRipple = max(abs(H_dB(passMask)));

                    if (actualStopMax <= -targetAs) && (actualPassRipple <= targetRp)
                        specMet = true;
                        app.AutoStatusLabel.Text = sprintf('Specs Met at Gen %d!', gen);
                        app.AutoStatusLabel.FontColor = [0.1 0.9 0.1]; % Green
                        break;
                    end
                end
            end

            % If user pressed stop manually
            if app.StopRequested && ~specMet
                app.AutoStatusLabel.Text = sprintf('Manually Stopped at Gen %d', gen);
                app.AutoStatusLabel.FontColor = [0.2 0.7 1.0]; % Light Blue
            end

            % Update Panel 3 Transfer Function Box
            if ~isempty(b_curr)
                tfString = app.formatTransferFunction(b_curr, a_curr);
                app.TFTextArea_2.Value = tfString;
            end

            % Reset Button States
            app.GenerateFilterButton.Enable = 'on';
            app.AutoSolveButton.Enable = 'on';
            app.StopAutoButton.Enable = 'off';
            app.IsOptimizing = false;
            app.StopRequested = false;
        
        
        
        end

        % Value changed function: FilterSelectDropDown_2
        function FilterSelectDropDown_2ValueChanged(app, event)
            updateSliderEdges_2(app);
            
        end

        % Value changed function: PassbandSlider_2
        function PassbandSlider_2ValueChanged(app, event)
            updateSliderEdges_2(app);
            
        end

        % Value changed function: StopBandSlider_2
        function StopBandSlider_2ValueChanged(app, event)
            updateSliderEdges_2(app);
            
        end

        % Button pushed function: StopAutoButton
        function StopAutoButtonPushed(app, event)
            if app.IsOptimizing
                app.StopRequested = true;
                app.AutoStatusLabel.Text = 'Status: Stopping gracefully...';
                app.AutoStatusLabel.FontColor = [1.0 0.4 0.2]; % Orange
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1145 735];
            app.UIFigure.Name = 'MATLAB App';

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.Title = 'Panel';
            app.Panel.Position = [1 0 465 736];

            % Create FilterSelectDropDownLabel
            app.FilterSelectDropDownLabel = uilabel(app.Panel);
            app.FilterSelectDropDownLabel.HorizontalAlignment = 'right';
            app.FilterSelectDropDownLabel.Position = [143 671 68 22];
            app.FilterSelectDropDownLabel.Text = 'Filter Select';

            % Create FilterSelectDropDown
            app.FilterSelectDropDown = uidropdown(app.Panel);
            app.FilterSelectDropDown.Items = {'LOW PASS', 'HIGH PASS', 'BAND PASS', 'BAND STOP'};
            app.FilterSelectDropDown.ValueChangedFcn = createCallbackFcn(app, @FilterSelectDropDownValueChanged, true);
            app.FilterSelectDropDown.Position = [226 671 98 22];
            app.FilterSelectDropDown.Value = 'LOW PASS';

            % Create FIRORIIRDropDownLabel
            app.FIRORIIRDropDownLabel = uilabel(app.Panel);
            app.FIRORIIRDropDownLabel.HorizontalAlignment = 'right';
            app.FIRORIIRDropDownLabel.Position = [143 635 71 22];
            app.FIRORIIRDropDownLabel.Text = 'FIR OR IIR?';

            % Create FIRORIIRDropDown
            app.FIRORIIRDropDown = uidropdown(app.Panel);
            app.FIRORIIRDropDown.Items = {'FIR', 'IIR'};
            app.FIRORIIRDropDown.Position = [229 635 56 22];
            app.FIRORIIRDropDown.Value = 'FIR';

            % Create GAORPSODropDownLabel
            app.GAORPSODropDownLabel = uilabel(app.Panel);
            app.GAORPSODropDownLabel.HorizontalAlignment = 'right';
            app.GAORPSODropDownLabel.Position = [143 599 78 22];
            app.GAORPSODropDownLabel.Text = 'GA OR PSO?';

            % Create GAORPSODropDown
            app.GAORPSODropDown = uidropdown(app.Panel);
            app.GAORPSODropDown.Items = {'GA', 'PSO'};
            app.GAORPSODropDown.Position = [236 599 59 22];
            app.GAORPSODropDown.Value = 'GA';

            % Create SIRORDERSpinnerLabel
            app.SIRORDERSpinnerLabel = uilabel(app.Panel);
            app.SIRORDERSpinnerLabel.HorizontalAlignment = 'right';
            app.SIRORDERSpinnerLabel.Position = [143 564 78 22];
            app.SIRORDERSpinnerLabel.Text = 'SIR ORDER?';

            % Create SIRORDERSpinner
            app.SIRORDERSpinner = uispinner(app.Panel);
            app.SIRORDERSpinner.Limits = [1 30];
            app.SIRORDERSpinner.RoundFractionalValues = 'on';
            app.SIRORDERSpinner.HorizontalAlignment = 'left';
            app.SIRORDERSpinner.Position = [236 564 95 22];
            app.SIRORDERSpinner.Value = 4;

            % Create PassbandSliderLabel
            app.PassbandSliderLabel = uilabel(app.Panel);
            app.PassbandSliderLabel.HorizontalAlignment = 'right';
            app.PassbandSliderLabel.Position = [42 477 62 22];
            app.PassbandSliderLabel.Text = 'Passband ';

            % Create PassbandSlider
            app.PassbandSlider = uislider(app.Panel, 'range');
            app.PassbandSlider.Limits = [0 1];
            app.PassbandSlider.ValueChangedFcn = createCallbackFcn(app, @PassbandSliderValueChanged, true);
            app.PassbandSlider.Position = [126 486 214 3];
            app.PassbandSlider.Value = [0 1];

            % Create StopBandSliderLabel
            app.StopBandSliderLabel = uilabel(app.Panel);
            app.StopBandSliderLabel.HorizontalAlignment = 'right';
            app.StopBandSliderLabel.Position = [43 425 61 22];
            app.StopBandSliderLabel.Text = 'Stop Band';

            % Create StopBandSlider
            app.StopBandSlider = uislider(app.Panel, 'range');
            app.StopBandSlider.Limits = [0 1];
            app.StopBandSlider.ValueChangedFcn = createCallbackFcn(app, @StopBandSliderValueChanged, true);
            app.StopBandSlider.Position = [126 434 214 3];
            app.StopBandSlider.Value = [0 1];

            % Create GenerateFilterButton
            app.GenerateFilterButton = uibutton(app.Panel, 'push');
            app.GenerateFilterButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateFilterButtonPushed, true);
            app.GenerateFilterButton.BackgroundColor = [0 1 0];
            app.GenerateFilterButton.FontName = 'Bookman Old Style';
            app.GenerateFilterButton.FontWeight = 'bold';
            app.GenerateFilterButton.FontAngle = 'italic';
            app.GenerateFilterButton.FontColor = [0 0 0];
            app.GenerateFilterButton.Position = [154 350 147 42];
            app.GenerateFilterButton.Text = 'Generate Filter';

            % Create Generation050Label
            app.Generation050Label = uilabel(app.Panel);
            app.Generation050Label.HorizontalAlignment = 'center';
            app.Generation050Label.Position = [147 312 162 22];
            app.Generation050Label.Text = 'Generation: 0 / 50';

            % Create StabilityNALabel
            app.StabilityNALabel = uilabel(app.Panel);
            app.StabilityNALabel.HorizontalAlignment = 'center';
            app.StabilityNALabel.Position = [137 288 187 22];
            app.StabilityNALabel.Text = 'Stability: N/A';

            % Create TFTextAreaLabel
            app.TFTextAreaLabel = uilabel(app.Panel);
            app.TFTextAreaLabel.HorizontalAlignment = 'right';
            app.TFTextAreaLabel.FontSize = 15;
            app.TFTextAreaLabel.Position = [17 222 63 37];
            app.TFTextAreaLabel.Text = {'Transfer'; 'Function'};

            % Create TFTextArea
            app.TFTextArea = uitextarea(app.Panel);
            app.TFTextArea.FontName = 'Monospaced';
            app.TFTextArea.FontSize = 15;
            app.TFTextArea.Position = [96 208 358 55];

            % Create MaxGenSpinnerLabel
            app.MaxGenSpinnerLabel = uilabel(app.Panel);
            app.MaxGenSpinnerLabel.HorizontalAlignment = 'right';
            app.MaxGenSpinnerLabel.Position = [143 530 67 22];
            app.MaxGenSpinnerLabel.Text = 'MAX GEN?';

            % Create MaxGenSpinner
            app.MaxGenSpinner = uispinner(app.Panel);
            app.MaxGenSpinner.RoundFractionalValues = 'on';
            app.MaxGenSpinner.HorizontalAlignment = 'left';
            app.MaxGenSpinner.Position = [225 530 122 22];
            app.MaxGenSpinner.Value = 50;

            % Create Panel2
            app.Panel2 = uipanel(app.UIFigure);
            app.Panel2.Title = 'Panel2';
            app.Panel2.Position = [833 1 313 735];

            % Create UIAxes
            app.UIAxes = uiaxes(app.Panel2);
            title(app.UIAxes, 'Magnitude Frequency Response (dB)')
            xlabel(app.UIAxes, 'Normalized Frequency (\times \pi rad/sample)')
            ylabel(app.UIAxes, 'Magnitude |H(e^{j\omega})|')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.YLim = [-80 5];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [1 515 300 185];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.Panel2);
            title(app.UIAxes2, 'Pole-Zero Map (Z-Plane)')
            xlabel(app.UIAxes2, 'Real Part')
            ylabel(app.UIAxes2, 'Imaginary Part')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.XLim = [-1.5 1.5];
            app.UIAxes2.YLim = [-1.5 1.5];
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Position = [1 235 300 185];

            % Create Panel3
            app.Panel3 = uipanel(app.UIFigure);
            app.Panel3.Title = 'Panel3';
            app.Panel3.Position = [465 -1 369 738];

            % Create FilterSelectDropDown_2Label
            app.FilterSelectDropDown_2Label = uilabel(app.Panel3);
            app.FilterSelectDropDown_2Label.HorizontalAlignment = 'right';
            app.FilterSelectDropDown_2Label.Position = [72 673 68 22];
            app.FilterSelectDropDown_2Label.Text = 'Filter Select';

            % Create FilterSelectDropDown_2
            app.FilterSelectDropDown_2 = uidropdown(app.Panel3);
            app.FilterSelectDropDown_2.Items = {'LOW PASS', 'HIGH PASS', 'BAND PASS', 'BAND STOP'};
            app.FilterSelectDropDown_2.ValueChangedFcn = createCallbackFcn(app, @FilterSelectDropDown_2ValueChanged, true);
            app.FilterSelectDropDown_2.Position = [155 673 98 22];
            app.FilterSelectDropDown_2.Value = 'LOW PASS';

            % Create FIRORIIRDropDown_2Label
            app.FIRORIIRDropDown_2Label = uilabel(app.Panel3);
            app.FIRORIIRDropDown_2Label.HorizontalAlignment = 'right';
            app.FIRORIIRDropDown_2Label.Position = [72 637 71 22];
            app.FIRORIIRDropDown_2Label.Text = 'FIR OR IIR?';

            % Create FIRORIIRDropDown_2
            app.FIRORIIRDropDown_2 = uidropdown(app.Panel3);
            app.FIRORIIRDropDown_2.Items = {'FIR', 'IIR'};
            app.FIRORIIRDropDown_2.Position = [158 637 56 22];
            app.FIRORIIRDropDown_2.Value = 'FIR';

            % Create GAORPSODropDown_2Label
            app.GAORPSODropDown_2Label = uilabel(app.Panel3);
            app.GAORPSODropDown_2Label.HorizontalAlignment = 'right';
            app.GAORPSODropDown_2Label.Position = [72 601 78 22];
            app.GAORPSODropDown_2Label.Text = 'GA OR PSO?';

            % Create GAORPSODropDown_2
            app.GAORPSODropDown_2 = uidropdown(app.Panel3);
            app.GAORPSODropDown_2.Items = {'GA', 'PSO'};
            app.GAORPSODropDown_2.Position = [165 601 59 22];
            app.GAORPSODropDown_2.Value = 'GA';

            % Create StopbandAttenEditFieldLabel
            app.StopbandAttenEditFieldLabel = uilabel(app.Panel3);
            app.StopbandAttenEditFieldLabel.HorizontalAlignment = 'right';
            app.StopbandAttenEditFieldLabel.Position = [30 566 136 22];
            app.StopbandAttenEditFieldLabel.Text = 'Min Stopband Atten (dB)';

            % Create StopbandAttenEditField
            app.StopbandAttenEditField = uieditfield(app.Panel3, 'numeric');
            app.StopbandAttenEditField.Position = [181 566 100 22];
            app.StopbandAttenEditField.Value = 35;

            % Create PassbandRippleEditFieldLabel
            app.PassbandRippleEditFieldLabel = uilabel(app.Panel3);
            app.PassbandRippleEditFieldLabel.HorizontalAlignment = 'right';
            app.PassbandRippleEditFieldLabel.Position = [25 533 148 22];
            app.PassbandRippleEditFieldLabel.Text = 'Max Passband Ripple (dB)';

            % Create PassbandRippleEditField
            app.PassbandRippleEditField = uieditfield(app.Panel3, 'numeric');
            app.PassbandRippleEditField.Position = [188 533 100 22];
            app.PassbandRippleEditField.Value = 0.5;

            % Create AutoSolveButton
            app.AutoSolveButton = uibutton(app.Panel3, 'push');
            app.AutoSolveButton.ButtonPushedFcn = createCallbackFcn(app, @AutoSolveButtonPushed, true);
            app.AutoSolveButton.BackgroundColor = [0 1 0];
            app.AutoSolveButton.FontName = 'Bookman Old Style';
            app.AutoSolveButton.FontWeight = 'bold';
            app.AutoSolveButton.FontAngle = 'italic';
            app.AutoSolveButton.FontColor = [0 0 0];
            app.AutoSolveButton.Position = [18 272 147 42];
            app.AutoSolveButton.Text = 'Generate Filter';

            % Create CurrentOrderLabel
            app.CurrentOrderLabel = uilabel(app.Panel3);
            app.CurrentOrderLabel.Position = [125 370 157 22];
            app.CurrentOrderLabel.Text = {'Estimated Order: --'; ''};

            % Create AutoStatusLabel
            app.AutoStatusLabel = uilabel(app.Panel3);
            app.AutoStatusLabel.Position = [60 334 230 22];
            app.AutoStatusLabel.Text = 'Status: Ready for Auto-Solve';

            % Create PassbandSlider_2Label
            app.PassbandSlider_2Label = uilabel(app.Panel3);
            app.PassbandSlider_2Label.HorizontalAlignment = 'right';
            app.PassbandSlider_2Label.Position = [25 489 62 22];
            app.PassbandSlider_2Label.Text = 'Passband ';

            % Create PassbandSlider_2
            app.PassbandSlider_2 = uislider(app.Panel3, 'range');
            app.PassbandSlider_2.Limits = [0 1];
            app.PassbandSlider_2.ValueChangedFcn = createCallbackFcn(app, @PassbandSlider_2ValueChanged, true);
            app.PassbandSlider_2.Position = [109 498 214 3];
            app.PassbandSlider_2.Value = [0 1];

            % Create StopBandSlider_2Label
            app.StopBandSlider_2Label = uilabel(app.Panel3);
            app.StopBandSlider_2Label.HorizontalAlignment = 'right';
            app.StopBandSlider_2Label.Position = [22 426 61 22];
            app.StopBandSlider_2Label.Text = 'Stop Band';

            % Create StopBandSlider_2
            app.StopBandSlider_2 = uislider(app.Panel3, 'range');
            app.StopBandSlider_2.Limits = [0 1];
            app.StopBandSlider_2.ValueChangedFcn = createCallbackFcn(app, @StopBandSlider_2ValueChanged, true);
            app.StopBandSlider_2.Position = [105 435 214 3];
            app.StopBandSlider_2.Value = [0 1];

            % Create TFTextAreaLabel_2
            app.TFTextAreaLabel_2 = uilabel(app.Panel3);
            app.TFTextAreaLabel_2.HorizontalAlignment = 'right';
            app.TFTextAreaLabel_2.FontSize = 15;
            app.TFTextAreaLabel_2.Position = [8 145 63 37];
            app.TFTextAreaLabel_2.Text = {'Transfer'; 'Function'};

            % Create TFTextArea_2
            app.TFTextArea_2 = uitextarea(app.Panel3);
            app.TFTextArea_2.FontName = 'Monospaced';
            app.TFTextArea_2.FontSize = 15;
            app.TFTextArea_2.Position = [87 131 266 55];

            % Create Generation0Label
            app.Generation0Label = uilabel(app.Panel3);
            app.Generation0Label.HorizontalAlignment = 'center';
            app.Generation0Label.Position = [188 292 162 22];
            app.Generation0Label.Text = 'Generation: 0 ';

            % Create StabilityNALabel_2
            app.StabilityNALabel_2 = uilabel(app.Panel3);
            app.StabilityNALabel_2.HorizontalAlignment = 'center';
            app.StabilityNALabel_2.Position = [176 268 187 22];
            app.StabilityNALabel_2.Text = 'Stability: N/A';

            % Create StopAutoButton
            app.StopAutoButton = uibutton(app.Panel3, 'push');
            app.StopAutoButton.ButtonPushedFcn = createCallbackFcn(app, @StopAutoButtonPushed, true);
            app.StopAutoButton.BackgroundColor = [0.851 0.2 0.2];
            app.StopAutoButton.FontName = 'Bookman Old Style';
            app.StopAutoButton.FontWeight = 'bold';
            app.StopAutoButton.FontAngle = 'italic';
            app.StopAutoButton.FontColor = [0 0 0];
            app.StopAutoButton.Enable = 'off';
            app.StopAutoButton.Position = [101 209 147 42];
            app.StopAutoButton.Text = 'STOP';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = dsp_filter_designer

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
