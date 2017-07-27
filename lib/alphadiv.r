
plot.alphadiv <- function(map00, alpha, metric) 
{
    p <- plot.boxplot.by.group.x.bmi(map00, alpha[metric], metric)
    save_plot(paste0("boxplot.", metric, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}


