 

process_scRNA_data <- function(raw_data_path, output_dir = "output_01", 
                               fig_width = 10, fig_height = 8, mtpct=5,dpi = 300) {
  #### 导入必要的包 ####
  required_packages <- c("Seurat", "tidyverse", "scutilsR", "ggplot2", 
                         "gprofiler2", "qs", "pbapply")
  for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
      stop(paste0("请安装 ", pkg, " 包"))
    }
  }
  
  #### 创建输出目录 ####
  dir.create(output_dir, showWarnings = FALSE)
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE)
  
  #### 数据导入和合并 ####
  # 获取所有样本文件夹
  sample_dirs <- list.dirs(raw_data_path, full.names = TRUE, recursive = FALSE)
  
  # 读取并处理每个样本
  seu_list <- lapply(sample_dirs, function(dir) {
    # 获取样本名称
    sample_name <- basename(dir)
    
    # 读取10x数据
    data <- Read10X(data.dir = dir)
    
    # 创建Seurat对象
    seu <- CreateSeuratObject(
      counts = data,
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
    
    # 添加样本信息
    seu$orig.ident <- sample_name
    seu$dornor <- sample_name
    
    return(seu)
  })
  
  # 合并所有样本
  if(length(seu_list) > 1) {
    seu <- merge(seu_list[[1]], 
                 y = seu_list[-1], 
                 add.cell.ids = names(seu_list),
                 project = "combined")
  } else {
    seu <- seu_list[[1]]
  }
  
  # 保存原始数据
  qs::qsave(seu, file.path(output_dir, "endo_raw.qs"))
  
  #### 质控部分 ####
  # 计算线粒体比例
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  
  # 质控可视化并保存图片
  # 基因数可视化
  gene_plot <- VlnPlot(seu, features = c("nFeature_RNA"), 
                       group.by = "orig.ident", pt.size = 0, raster=FALSE) + 
    NoLegend() + 
    geom_hline(yintercept = c(300,7500), linetype = "dashed") +
    ggtitle("Gene Counts Distribution")
  
  ggsave(file.path(fig_dir, "gene_counts.pdf"), 
         gene_plot, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "gene_counts.png"), 
         gene_plot, width = fig_width, height = fig_height, dpi = dpi)
  
  # UMI数可视化
  umi_plot <- VlnPlot(seu, features = c("nCount_RNA"), 
                      group.by = "orig.ident", pt.size = 0, raster=FALSE, log = TRUE) + 
    NoLegend() + 
    geom_hline(yintercept = c(2000,40000), linetype = "dashed") +
    ggtitle("UMI Counts Distribution")
  
  ggsave(file.path(fig_dir, "umi_counts.pdf"), 
         umi_plot, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "umi_counts.png"), 
         umi_plot, width = fig_width, height = fig_height, dpi = dpi)
  
  # 线粒体比例可视化
  mt_plot <- VlnPlot(seu, features = c("percent.mt"), 
                     group.by = "orig.ident", pt.size = 0, raster=FALSE, log = TRUE) + 
    NoLegend() + 
    geom_hline(yintercept = c(mtpct), linetype = "dashed") +
    ggtitle("Mitochondrial Percentage")
  
  ggsave(file.path(fig_dir, "mt_percent.pdf"), 
         mt_plot, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "mt_percent.png"), 
         mt_plot, width = fig_width, height = fig_height, dpi = dpi)
  
  # 标记质控结果
  seu$QC <- "high"
  seu$QC[seu$nFeature_RNA < 300 | seu$percent.mt > 5 ] <- "low"
  
  # 质控散点图
  qc_scatter <- ggplot(seu@meta.data, aes(percent.mt, nFeature_RNA, color = QC)) +
    geom_point(size = 1) +
    geom_density_2d(color = "#364f6b", contour_var = "ndensity") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("#f08a5d", "#364f6b")) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    facet_wrap(~dornor, ncol = 2) +
    ggtitle("QC Metrics Scatter Plot")
  
  ggsave(file.path(fig_dir, "qc_scatter.pdf"), 
         qc_scatter, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "qc_scatter.png"), 
         qc_scatter, width = fig_width, height = fig_height, dpi = dpi)
  
  # 保存质控结果表格
  qc_stats <- table(seu$orig.ident, seu$QC)
  write.csv(qc_stats, file.path(output_dir, "qc_statistics.csv"))
  
  # 过滤低质量细胞
  seu <- subset(seu, QC == "high")
  qs::qsave(seu, file.path(output_dir, "endo_QC.qs"))
  
  #### 标准化和降维 ####
  # 标准化
  seu <- NormalizeData(seu)
  
  # 找到高变基因
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  
  # 数据缩放
  seu <- ScaleData(seu, features = rownames(seu))
  
  # PCA分析
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  
  # UMAP降维
  seu <- RunUMAP(seu, dims = 1:30)
  
  #### 细胞周期评分 ####
  # 获取小鼠同源基因
  mmus_s = gorth(cc.genes.updated.2019$s.genes, 
                 source_organism = "hsapiens", 
                 target_organism = "mmusculus")$ortholog_name
  mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, 
                   source_organism = "hsapiens", 
                   target_organism = "mmusculus")$ortholog_name
  
  seu <- CellCycleScoring(seu, s.features = mmus_s, g2m.features = mmus_g2m)
  
  # 细胞周期可视化
  cell_cycle_plot <- DimPlot(seu, group.by = "Phase") +
    ggtitle("Cell Cycle Phase Distribution")
  
  ggsave(file.path(fig_dir, "cell_cycle.pdf"), 
         cell_cycle_plot, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "cell_cycle.png"), 
         cell_cycle_plot, width = fig_width, height = fig_height, dpi = dpi)
  
  qs::qsave(seu, file.path(output_dir, "endo_QC.qsCC.qs"))
  
  #### 双细胞标注 ####
  seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident", PCs = 1:10)
  
  # 双细胞标注可视化
  doublet_plot <- DimPlot(seu, group.by = "DF.classifications") +
    ggtitle("Doublet Classifications")
  
  ggsave(file.path(fig_dir, "doublets.pdf"), 
         doublet_plot, width = fig_width, height = fig_height, dpi = dpi)
  ggsave(file.path(fig_dir, "doublets.png"), 
         doublet_plot, width = fig_width, height = fig_height, dpi = dpi)
  
  # 保存双细胞统计结果
  doublet_stats <- table(seu$DF.classifications)
  write.csv(as.data.frame(doublet_stats), 
            file.path(output_dir, "doublet_statistics.csv"))
  
  qs::qsave(seu, file.path(output_dir, "endo_QC_CC_markdup.qs"))
  
  # 返回处理后的对象和结果路径
  return(list(
    seurat_object = seu,
    output_directory = output_dir,
    figure_directory = fig_dir,
    qc_statistics = qc_stats,
    doublet_statistics = doublet_stats
  ))
}

# 运行函数
results <- process_scRNA_data(
  raw_data_path = "/mnt/i/scRNA-seq_WHJ_SYC/count",
  output_dir = "/mnt/i/scRNA-seq_WHJ_SYC/output_01",
  fig_width = 12,    # 可以调整图片宽度 
  fig_height = 8,    # 可以调整图片高度
  dpi = 300         # 可以调整图片分辨率
)

# 查看结果目录
list.files(results$output_directory)
list.files(results$figure_directory)




#### 导入数据 ####
samples <- list.dirs("",full.names = F,recursive = F)
samples   ##### 目录是matrix文件目录  下面的基因名一定要指定第一列，默认是2
seu.list <- pbapply::pblapply(samples,function(sn) {
  counts <- Read10X(file.path("/mnt/g/desktop/endo/raw", sn))
  sn = gsub("_","-",sn)
  colnames(counts) <- paste(sn,colnames(counts),sep = "_")
  seu <- CreateSeuratObject(counts = counts,)
  return(seu)
}
)

seu.list

## 合并样本
seu <- base::Reduce(f = merge, x = seu.list)
seu
seu$dornor <- sub("\\..*", "", seu$orig.ident)
seu$group <- sub("-.*", "", seu$dornor)
## 清理不要的变量 
rm(seu.list)
gc()



Idents(art) <- "group"
seusss <- FindAllMarkers(seu,assay = "RNA",test.use = "MAST",logfc.threshold = 0,min.pct = 0,only.pos = T)

 