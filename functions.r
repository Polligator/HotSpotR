CreateHS <- function(object, model = "normal", tree = F, num_umi = NULL, logdata = FALSE) {
    library(reticulate)
    library(Seurat)
    hotspot <- import("hotspot", convert = F)
    warnings <- import("warnings")
    warnings$filterwarnings("ignore")
    workers <- getOption("mc.cores")
    
    if (is.null(workers)) {
        workers <- 1
    }
    exprData = log1p(object@assays$RNA@counts)

    exprData = as.data.frame(as.matrix(exprData))

    # remove genes that do not have any standard deviation
    sds = apply(exprData, 1, sd)
    exprData = exprData[which(sds > 0), ]

    # generate the Hotspot object in python, potentially using the tree
    if (tree) {
        message("Using Tree,there is no tree")
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot$legacy_init(exprData, tree = pyTree, model = model)
        } else {
            py$umi_df <- r_to_py(data.frame(num_umi))
            py_run_string("umi_counts = umi_df.iloc[:, 0]")
            hs <- hotspot$Hotspot$legacy_init(exprData, tree = pyTree, model = model, umi_counts = py$umi_counts)
        }
    } else {
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot$legacy_init(exprData, latent = as.data.frame(Embeddings(object, reduction = "pca")), model = model)
        } else {
            py$umi_df <- r_to_py(data.frame(num_umi))
            py_run_string("umi_counts = umi_df.iloc[:, 0]")
            hs <- hotspot$Hotspot$legacy_init(exprData, latent = as.data.frame(Embeddings(object, reduction = "pca")), model = model, umi_counts = py$umi_counts)
        }
    }
    return(hs)
}


CreateKnn <- function(hs,n_neighbors=NULL) {
  # create knn graph, specify nn or use object default
  if (is.null(n_neighbors)) {
    "missing n_neighbors value"
  } else {
    hs$create_knn_graph(F, n_neighbors = as.integer(n_neighbors))
  }
  return(hs)
}
