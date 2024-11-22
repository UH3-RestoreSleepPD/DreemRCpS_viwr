% for (i in 1:length(years_to_draw)) {
%     curr_year <- curr[curr$YEAR == years_to_draw[i],]
%     if (i == length(years_to_draw)) {
%         col <- "#b31dc2"
%         lwd <- max_lwd
%         } else {
%             col <- "#cd24de"
%             lwd <- max_lwd * 1 / ((length(years_to_draw)-i+.5)^1.5)
%             }
% 
%         lines(curr_year$AGE, 100*curr_year$prop, col=col, lwd=lwd)
%         }



