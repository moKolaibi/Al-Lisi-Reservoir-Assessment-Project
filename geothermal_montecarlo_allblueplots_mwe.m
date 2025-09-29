function geothermal_montecarlo_allblueplots_mwe(type, data, num_bins, stats, cis, sample_size)
% Publication-ready plotting for Monte Carlo geothermal simulation (MWe)
% Improvements: all annotations inside figure, most likely interval patch height matches histogram, large/bold fonts

if nargin < 6, sample_size = []; end

switch lower(type)
    case 'pdf'
        figure('Color','w','Name','PDF of Electrical Power Harvested','NumberTitle','off');
        h_hist = histogram(data, num_bins, 'Normalization', 'pdf', ...
            'FaceColor', [0.16 0.44 0.74], 'EdgeColor', [0.08 0.22 0.37], ...
            'LineWidth',2, 'FaceAlpha',0.85, 'DisplayName','Power Distribution');
        xlabel('Electrical Power Harvested (MWe)', 'FontSize',18,'FontWeight','bold');
        ylabel('Probability Density', 'FontSize',18,'FontWeight','bold');
        title({'PDF of Electrical Power Harvested', ...
            sprintf('Monte Carlo Simulation, N = %d', sample_size)}, 'FontSize',20,'FontWeight','bold');
        set(gca,'FontSize',15,'FontWeight','bold','Box','on','LineWidth',1.5);
        grid on;

        ylimv = ylim;
        xlimv = xlim;

        % Highlight mode bin (most likely interval) - patch matches bin height
        [~, mode_idx] = max(h_hist.Values);
        x_mode = [h_hist.BinEdges(mode_idx) h_hist.BinEdges(mode_idx+1)];
        mode_height = h_hist.Values(mode_idx);
        hold on;
        h_patch = patch([x_mode(1) x_mode(2) x_mode(2) x_mode(1)], ...
                        [0 0 mode_height mode_height], [0.85 0.54 0.25], ...
                        'FaceAlpha',0.7, 'EdgeColor','none', 'DisplayName','Most Likely Interval');

        % Most likely interval annotation (above mode bin, inside plot)
        y_mode = mode_height * 1.08;
        x_mode_center = mean(x_mode);
        text(x_mode_center, y_mode, ...
            sprintf('Most Likely Bin\n%.2f -- %.2f MWe\nMode: %.2f MWe', x_mode(1), x_mode(2), x_mode_center), ...
            'Color', [0.3 0.3 0.15], 'FontSize',16, 'FontWeight','bold', ...
            'HorizontalAlignment','center','BackgroundColor','w','Margin',3);

        % CI vertical lines and dots at base
        plot([cis(1) cis(1)], ylimv, 'r-', 'LineWidth',2, 'HandleVisibility','off');
        plot([cis(2) cis(2)], ylimv, 'r-', 'LineWidth',2, 'HandleVisibility','off');
        plot(cis(1), 0, 'ro', 'MarkerFaceColor','r','MarkerSize',10, 'HandleVisibility','off');
        plot(cis(2), 0, 'ro', 'MarkerFaceColor','r','MarkerSize',10, 'HandleVisibility','off');

        % CI labels INSIDE plot (slightly above x-axis)
        y_cilabel = ylimv(2)*0.08;
        text(cis(1), y_cilabel, sprintf('10%% CI: %.2f MWe', cis(1)), ...
            'Color','r','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w');
        text(cis(2), y_cilabel, sprintf('90%% CI: %.2f MWe', cis(2)), ...
            'Color','r','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w');

        % Stats box INSIDE plot (top right corner)
        x_stats = xlimv(1) + 0.78*diff(xlimv);
        y_stats = ylimv(2)*0.94;
        stats_str = {
            sprintf('Min:   %.2f MWe', stats(1))
            sprintf('Max:   %.2f MWe', stats(2))
            sprintf('Mode:  %.2f MWe', stats(3))
            sprintf('Std:   %.2f MWe', stats(4))
            sprintf('10%% CI: %.2f MWe', cis(1))
            sprintf('90%% CI: %.2f MWe', cis(2))
            };
        text(x_stats, y_stats, strjoin(stats_str,'\n'), ...
            'FontSize',14,'FontName','Consolas','Color',[0.1 0.3 0.5], ...
            'BackgroundColor','w','FontWeight','bold','Margin',4,'VerticalAlignment','top');

        legend([h_hist h_patch], {'Power Distribution','Most Likely Interval'},'Location','northeast','FontSize',15);
        hold off;

    case 'cdf'
        figure('Color','w','Name','CDF of Electrical Power Harvested','NumberTitle','off');
        % Create empirical CDF for exact percentiles
        [sorted_data, sort_idx] = sort(data);
        y_cdf = (1:length(sorted_data))/length(sorted_data);

        h_line = plot(sorted_data, y_cdf, 'LineWidth',3, 'Color', [0 0.3 0.8], 'DisplayName','Cumulative Probability');
        xlabel('Electrical Power Harvested (MWe)', 'FontSize',18,'FontWeight','bold');
        ylabel('Cumulative Probability', 'FontSize',18,'FontWeight','bold');
        title({'CDF of Electrical Power Harvested', ...
            sprintf('Monte Carlo Simulation, N = %d', sample_size)}, 'FontSize',20,'FontWeight','bold');
        set(gca,'FontSize',15,'FontWeight','bold','Box','on','LineWidth',1.5);
        grid on; hold on;

        xlimv = xlim; ylimv = ylim;

        % Percentile markers INSIDE plot -- use exact input stats values!
        colormap = {[0 0.7 0], [0 0.2 1], [0.8 0 0]}; % green, blue, red
        labels = {'P_{10}', 'P_{50}', 'P_{90}'};
        h_dots = gobjects(3,1);
        percent_values = [0.10, 0.50, 0.90];
        for k=1:3
            x_val = stats(k); % Use exact P10, P50, P90
            % Find corresponding y_val from empirical CDF
            y_val = interp1(sorted_data, y_cdf, x_val, 'linear','extrap');
            h_dots(k) = plot(x_val, y_val, 'o', 'MarkerSize',13, ...
                'MarkerFaceColor',colormap{k}, 'MarkerEdgeColor','k', 'LineWidth',2, 'DisplayName',labels{k});
            % Place annotation just above the marker
            text(x_val, y_val + 0.07, sprintf('%s: %.2f MWe', labels{k}, x_val), ...
                'Color',colormap{k}, 'FontSize',17, 'FontWeight','bold', ...
                'BackgroundColor','w','HorizontalAlignment','center','VerticalAlignment','bottom');
        end

        % Stats box INSIDE plot (top right corner)
        x_stats = xlimv(1) + 0.78*diff(xlimv);
        y_stats = ylimv(2)*0.94;
        stats_str = {
            sprintf('Min:   %.2f MWe', min(data))
            sprintf('Max:   %.2f MWe', max(data))
            sprintf('Mean:  %.2f MWe', mean(data))
            sprintf('Std:   %.2f MWe', std(data))
            sprintf('P_{10}: %.2f MWe', stats(1))
            sprintf('P_{50}: %.2f MWe', stats(2))
            sprintf('P_{90}: %.2f MWe', stats(3))
        };
        text(x_stats, y_stats, strjoin(stats_str,'\n'), ...
            'FontSize',14,'FontName','Consolas','Color',[0.1 0.3 0.5], ...
            'BackgroundColor','w','FontWeight','bold','Margin',4,'VerticalAlignment','top');

        legend([h_line h_dots'], {'Cumulative Probability','P_{10}','P_{50}','P_{90}'}, 'FontSize',16, ...
            'Location','northeast','Box','off');
        hold off;

    case 'thermal'
        figure('Color','w','Name','PDF of Total Thermal Energy','NumberTitle','off');
        h_hist = histogram(data, num_bins, 'Normalization', 'pdf', ...
            'FaceColor', [0.16 0.44 0.74], 'EdgeColor', [0.08 0.22 0.37], 'LineWidth',2,'FaceAlpha',0.80, ...
            'DisplayName','Thermal Energy Distribution');
        xlabel('Total Thermal Energy (J)', 'FontSize',18,'FontWeight','bold');
        ylabel('Probability Density', 'FontSize',18,'FontWeight','bold');
        title({'PDF of Total Thermal Energy', ...
            sprintf('Monte Carlo Simulation, N = %d', sample_size)}, 'FontSize',20,'FontWeight','bold');
        set(gca,'FontSize',15,'FontWeight','bold','Box','on','LineWidth',1.5);
        grid on;

        ylimv = ylim; xlimv = xlim;
        [~, mode_idx] = max(h_hist.Values);
        x_mode = [h_hist.BinEdges(mode_idx) h_hist.BinEdges(mode_idx+1)];
        mode_height = h_hist.Values(mode_idx);
        hold on;
        h_patch = patch([x_mode(1) x_mode(2) x_mode(2) x_mode(1)], ...
                        [0 0 mode_height mode_height], [0.85 0.54 0.25], ...
                        'FaceAlpha',0.7, 'EdgeColor','none', 'DisplayName','Most Likely Interval');

        % Most likely interval annotation (above mode bin, inside plot)
        y_mode = mode_height * 1.08;
        x_mode_center = mean(x_mode);
        text(x_mode_center, y_mode, ...
            sprintf('Most Likely Bin\n%.2e -- %.2e J\nMode: %.2e J', x_mode(1), x_mode(2), x_mode_center), ...
            'Color', [0.3 0.3 0.15], 'FontSize',16, 'FontWeight','bold', ...
            'HorizontalAlignment','center','BackgroundColor','w','Margin',3);

        % CI vertical lines and dots
        plot([cis(1) cis(1)], ylimv, 'r-', 'LineWidth',2, 'HandleVisibility','off');
        plot([cis(2) cis(2)], ylimv, 'r-', 'LineWidth',2, 'HandleVisibility','off');
        plot(cis(1), 0, 'ro', 'MarkerFaceColor','r','MarkerSize',10, 'HandleVisibility','off');
        plot(cis(2), 0, 'ro', 'MarkerFaceColor','r','MarkerSize',10, 'HandleVisibility','off');

        % CI labels INSIDE plot (slightly above x-axis)
        y_cilabel = ylimv(2)*0.08;
        text(cis(1), y_cilabel, sprintf('10%% CI: %.2e J', cis(1)), ...
            'Color','r','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w');
        text(cis(2), y_cilabel, sprintf('90%% CI: %.2e J', cis(2)), ...
            'Color','r','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w');

        % Stats box INSIDE plot (top right corner)
        x_stats = xlimv(1) + 0.78*diff(xlimv);
        y_stats = ylimv(2)*0.94;
        stats_str = {
            sprintf('Min:   %.2e J', stats(1))
            sprintf('Max:   %.2e J', stats(2))
            sprintf('Mode:  %.2e J', stats(3))
            sprintf('Mean:  %.2e J', stats(4))
            sprintf('Std:   %.2e J', stats(5))
            sprintf('10%% CI: %.2e J', cis(1))
            sprintf('90%% CI: %.2e J', cis(2))
            };
        text(x_stats, y_stats, strjoin(stats_str,'\n'), ...
            'FontSize',14,'FontName','Consolas','Color',[0.1 0.3 0.5], ...
            'BackgroundColor','w','FontWeight','bold','Margin',4,'VerticalAlignment','top');

        legend([h_hist h_patch], {'Thermal Energy Distribution','Most Likely Interval'},'Location','northeast','FontSize',15);
        hold off;
end

set(gcf, 'Position', [100 100 1100 600]); % Wider figure for publication
end