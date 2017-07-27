get.pch <- function(map)
{    
    lookup <- c(19,17) # point type
    names(lookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch <- lookup[as.character(map$Ethnicity)] 
    return(pch)
}