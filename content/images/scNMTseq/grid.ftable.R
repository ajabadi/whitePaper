grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
    
    nc <- ncol(d)
    nr <- nrow(d)
    
    ## character table with added row and column names
    extended_matrix <- cbind(c("", rownames(d)),
                             rbind(colnames(d),
                                   as.matrix(d)))
    
    ## string width and height
    w <- apply(extended_matrix, 2, strwidth, "inch")
    h <- apply(extended_matrix, 2, strheight, "inch")
    
    widths <- apply(w, 2, max)
    heights <- apply(h, 1, max)*4
    
    padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)
    
    x <- cumsum(widths + padding) - 0.5 * padding
    y <- cumsum(heights + padding) - padding
    
    rg <- rectGrob(x = unit(x - widths/2, "in"),
                   y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                   width = unit(widths + padding, "in"),
                   height = unit(heights + padding, "in"))
    
    tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                   y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                   just = "center")
    
    g <- gTree(children = gList(rg, tg), ...,
               x = x, y = y, widths = widths, heights = heights)
    
    grid.draw(g)
}
