#' Planar Maximally Filtered Graph
#' @description Applies the Planar Maximally Filtered Graph (PMFG) filtering method
#' (see and cite Tumminello et al., 2005).
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous?
#' Defaults to FALSE.
#' Set to TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @param na.data How should missing data be handled? For "listwise" deletion the \emph{na.omit} fucntion is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a PMFG-filtered associaton matrix
#' @references
#' Carey V, Long L and Gentleman R (2017).
#' RBGL: An interface to the BOOST graph library.
#' R package version 1.54.0, http://www.bioconductor.org.
#'  
#' Tumminello, M., Aste, T., Di Matteo, T., & Mantegna, R. N. (2005).
#' A tool for filtering information in complex systems.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{102}(30), 10421-10426.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom utils installed.packages
#' @export
#PMFG Filtering Method----
PMFG <- function (data, binary = FALSE, weighted = TRUE,
                  na.data = c("listwise","fiml"), progBar = TRUE)
{
    
    if(!"RBGL" %in% rownames(installed.packages())){
        cat("In order to perform this function, please copy code below to install: RBGL and graph packages",sep="\n")
        cat('install.packages("BiocManager")',sep="\n")
        cat('BiocManager::install(version = "3.14")',sep="\n")
        cat('BiocManager::install("RBGL")',sep="\n")
    }
    
    #missing data handling
    if(missing(na.data)){
        if(any(is.na(data)))
        {stop("Missing values were detected! Set 'na.data' argument")
        }else{na.data<-"none"}
    }else{na.data<-match.arg(na.data)}
    
    if(na.data=="listwise"){
        rem<-na.action(na.omit(data))
        warning(paste(length(na.action(na.omit(data)))),
                " rows were removed for missing data\nrow(s): ",
                paste(na.action(na.omit(data)),collapse = ", "))
        data<-na.omit(data)
    }else if(na.data=="fiml"){
        data<-psych::corFiml(data)
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(binary){cormat<-psych::tetrachoric(data)$rho
    }else{cormat<-cor(data)}
    
    #create sparse data
    i<-as.vector(rep(1:ncol(data),ncol(data)))
    j<-sort(as.vector(rep(1:ncol(data),ncol(data))))
    w<-as.vector(cormat)
    
    kk<-which(i<j)
    
    ijw<-cbind(i[kk],j[kk],w[kk])
    
    ijw<-ijw[order(ijw[,3],decreasing=TRUE),]
    
    P<-Matrix::Matrix(0,nrow=ncol(data),ncol=ncol(data))
    
    as_graphnel <- function(graph) {
        
        if (!igraph::is_igraph(graph)) {
            stop("Not an igraph graph")
        }
        
        if ("name" %in% suppressWarnings(igraph::vertex_attr_names(graph)) &&
            is.character(suppressWarnings(igraph::V(graph)$name))) {
            name <- suppressWarnings(igraph::V(graph)$name)
        } else {
            name <- as.character(seq(igraph::vcount(graph)))    
        }
        
        edgemode <- "undirected"  
        
        if ("weight" %in% igraph::edge_attr_names(graph) &&
            is.numeric(igraph::E(graph)$weight)) {
            al <- lapply(igraph::as_adj_edge_list(graph, "out"), as.vector)
            for (i in seq(along=al)) {
                edges <- igraph::ends(graph, al[[i]], names = FALSE)
                edges <- ifelse( edges[,2]==i, edges[,1], edges[,2])
                weights <- igraph::E(graph)$weight[al[[i]]]
                al[[i]] <- list(edges=edges, weights=weights)
            }
        } else {
            al <- igraph::as_adj_list(graph, "out")
            al <- lapply(al, function(x) list(edges=as.vector(x)))
        }  
        
        names(al) <- name
        res <- graph::graphNEL(nodes=name, edgeL=al, edgemode=edgemode)
        
        ## Add graph attributes (other than 'directed')
        ## Are this "officially" supported at all?
        
        g.n <- igraph::graph_attr_names(graph)
        if ("directed" %in% g.n) {
            warning("Cannot add graph attribute `directed'")
            g.n <- g.n[ g.n != "directed" ]
        }
        for (n in g.n) {
            res@graphData[[n]] <- igraph::graph_attr(graph, n)
        }
        
        ## Add vertex attributes (other than 'name', that is already
        ## added as vertex names)
        
        v.n <- igraph::vertex_attr_names(graph)
        v.n <- v.n[ v.n != "name" ]
        for (n in v.n) {
            graph::nodeDataDefaults(res, attr=n) <- NA
            graph::nodeData(res, attr=n) <- igraph::vertex_attr(graph, n)
        }
        
        ## Add edge attributes (other than 'weight')
        
        e.n <- igraph::edge_attr_names(graph)
        e.n <- e.n[ e.n != "weight" ]
        if (length(e.n) > 0) {
            el <- igraph::as_edgelist(graph)
            el <- paste(sep="|", el[,1], el[,2])
            for (n in e.n) {
                graph::edgeDataDefaults(res, attr=n) <- NA
                res@edgeData@data[el] <- mapply(function(x,y) {
                    xx <- c(x,y); names(xx)[length(xx)] <- n; xx },
                    res@edgeData@data[el],
                    igraph::edge_attr(graph, n),
                    SIMPLIFY=FALSE)
            }
        }
        
        res
    }
    
    for(ii in 1:pmin(6,nrow(ijw))){
        P[ijw[ii,1],ijw[ii,2]]<-ijw[ii,3]
        P[ijw[ii,2],ijw[ii,1]]<-ijw[ii,3]
    }
    
    E<-6
    P1<-P
    
    if(progBar==TRUE){pb <- txtProgressBar(max=(3*(ncol(data)-2)), style = 3)}
    
    while(E < 3*(ncol(data)-2)){
        ii<-ii+1
        P1[ijw[ii,1],ijw[ii,2]]<-ijw[ii,3]
        P1[ijw[ii,2],ijw[ii,1]]<-ijw[ii,3]
        
        graph<-suppressWarnings(
            igraph::as.igraph(qgraph::qgraph(P1,DoNotPlot=TRUE))   
        )
        
        g<-as_graphnel(graph)
        
        if(RBGL::boyerMyrvoldPlanarityTest(g)==TRUE){
            P<-P1
            E<-E+1
            
            if(progBar==TRUE)
            {setTxtProgressBar(pb, E)}
        }else{P1<-P}
        
        if(ii>(ncol(data)*(ncol(data)-1)/2)){message("PMFG not found")}
        
    }
    if(progBar==TRUE)
    {close(pb)}
    
    pmfg<-as.matrix(P)
    
    if(!weighted)
    {pmfg<-ifelse(pmfg!=0,1,0)}
    
    return(pmfg)
}
#----