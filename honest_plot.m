function  [ med, quant, whisk, clean_data ] = honest_plot(data, varargin)
    opt	= ParseArgs(varargin,...
            'whiskers'		, 1.5);
   
    quant = quantile(data,[.25 .75]);
    med = median(data);
    crit(1,:) = quant(1,:)-(quant(2,:)-quant(1,:))*opt.whiskers;
    crit(2,:) = quant(2,:)+(quant(2,:)-quant(1,:))*opt.whiskers;
    outl = logical((data < crit(1,:)) +  (data > crit(2,:)));
    clean_data = data;
    out_data = data;
    clean_data(outl) = NaN;
    out_data(~outl) = NaN;
    % use most extreme data point as whisker value
    whisk(1,:) = nanmin(clean_data);
    whisk(2,:) = nanmax(clean_data);
    
    % do plotting
    x_vals = 1:size(data,2);
    hold on
    patch_x = [x_vals, fliplr(x_vals)];
    patch_y = [quant(1,:),fliplr(quant(2,:))]
    patch(patch_x,patch_y,'b', 'edgecolor', 'none');
    plot(x_vals, whisk,'b--');
    plot(x_vals, med,'wo-')
    plot(x_vals, out_data,'bx')
    
end

