# R脚本：基于基因相关性分析构建蛋白质相互作用网络
# 作者：AI助手
# 日期：2025-07-22

# 加载必要的包
library(data.table)
library(dplyr)
library(igraph)

# 设置文件路径
correlation_file <- "D:/biological_relevence_for_FM_feature/tryfile_FTV_T0_corr_res_spear_train.csv"
converted_links_file <- "D:/biological_relevence_for_FM_feature/9606.protein.links.v12.0_symbol_converted.txt"

# 读取相关性文件
cat("正在读取基因相关性文件...\n")
gene_correlation <- fread(correlation_file, header = TRUE)

# 检查相关性文件结构
cat("基因相关性文件结构：\n")
print(head(gene_correlation))
cat("相关性文件维度：", dim(gene_correlation), "\n")
cat("相关性文件列名：", names(gene_correlation), "\n")

# 标准化列名
if(ncol(gene_correlation) >= 2) {
  colnames(gene_correlation)[1:2] <- c("Gene_Name", "Correlation_Coefficient")
} else {
  stop("相关性文件应该至少包含两列：基因名和相关性系数")
}

# 读取转换后的protein连接文件
cat("正在读取转换后的protein连接文件...\n")
protein_links <- fread(converted_links_file, header = TRUE)

# 检查连接文件结构
cat("Protein连接文件结构：\n")
print(head(protein_links))
cat("连接文件维度：", dim(protein_links), "\n")
cat("连接文件列名：", names(protein_links), "\n")

# 标准化连接文件列名
protein1_col <- names(protein_links)[1]
protein2_col <- names(protein_links)[2]
if(ncol(protein_links) >= 3) {
  score_col <- names(protein_links)[3]
} else {
  score_col <- NULL
}

cat("使用的列名：", protein1_col, "，", protein2_col, "\n")
if(!is.null(score_col)) {
  cat("分数列：", score_col, "\n")
}

# 获取相关性文件中的基因列表
correlation_genes <- unique(gene_correlation$Gene_Name)
cat("相关性文件中的基因数量：", length(correlation_genes), "\n")

# 显示基因示例
cat("基因名示例：\n")
print(head(correlation_genes, 10))

# 获取protein连接文件中的所有基因
protein_genes <- unique(c(protein_links[[protein1_col]], protein_links[[protein2_col]]))
# 移除NA值
protein_genes <- protein_genes[!is.na(protein_genes)]
cat("Protein连接文件中的基因数量：", length(protein_genes), "\n")

# 找到两个文件中共同的基因
common_genes <- intersect(correlation_genes, protein_genes)
cat("共同基因数量：", length(common_genes), "\n")
cat("覆盖率（相关性基因在连接文件中的比例）：", 
    round(length(common_genes)/length(correlation_genes)*100, 2), "%\n")

# 显示共同基因示例
if(length(common_genes) > 0) {
  cat("共同基因示例：\n")
  print(head(common_genes, 10))
} else {
  cat("警告：没有找到共同的基因！请检查基因名格式是否一致。\n")
  cat("相关性文件基因示例：", head(correlation_genes, 5), "\n")
  cat("连接文件基因示例：", head(protein_genes, 5), "\n")
}

# 筛选出涉及相关性基因的protein连接
cat("正在筛选涉及相关性基因的protein连接...\n")
relevant_links <- protein_links[
  (protein_links[[protein1_col]] %in% correlation_genes) | 
  (protein_links[[protein2_col]] %in% correlation_genes)
]

cat("筛选后的连接数量：", nrow(relevant_links), "\n")
cat("原始连接数量：", nrow(protein_links), "\n")
cat("筛选比例：", round(nrow(relevant_links)/nrow(protein_links)*100, 2), "%\n")

# 进一步筛选：只保留两个基因都在相关性列表中的连接
cat("正在筛选两个基因都在相关性列表中的连接...\n")
both_genes_relevant <- relevant_links[
  (relevant_links[[protein1_col]] %in% correlation_genes) & 
  (relevant_links[[protein2_col]] %in% correlation_genes)
]

cat("两个基因都相关的连接数量：", nrow(both_genes_relevant), "\n")

# 为连接添加相关性信息
cat("正在为连接添加相关性信息...\n")

# 创建基因到相关性的映射
gene_to_corr <- setNames(gene_correlation$Correlation_Coefficient, gene_correlation$Gene_Name)

# 为relevant_links添加相关性信息
relevant_links_with_corr <- relevant_links
relevant_links_with_corr$Gene1_Correlation <- gene_to_corr[relevant_links_with_corr[[protein1_col]]]
relevant_links_with_corr$Gene2_Correlation <- gene_to_corr[relevant_links_with_corr[[protein2_col]]]

# 为both_genes_relevant添加相关性信息
both_genes_with_corr <- both_genes_relevant
both_genes_with_corr$Gene1_Correlation <- gene_to_corr[both_genes_with_corr[[protein1_col]]]
both_genes_with_corr$Gene2_Correlation <- gene_to_corr[both_genes_with_corr[[protein2_col]]]

# 计算相关性统计
if(nrow(both_genes_with_corr) > 0) {
  cat("相关性统计（两个基因都相关的连接）：\n")
  cat("基因1相关性范围：", range(both_genes_with_corr$Gene1_Correlation, na.rm = TRUE), "\n")
  cat("基因2相关性范围：", range(both_genes_with_corr$Gene2_Correlation, na.rm = TRUE), "\n")
  cat("平均相关性（基因1）：", round(mean(both_genes_with_corr$Gene1_Correlation, na.rm = TRUE), 4), "\n")
  cat("平均相关性（基因2）：", round(mean(both_genes_with_corr$Gene2_Correlation, na.rm = TRUE), 4), "\n")
}

# 保存结果文件
output_dir <- "D:/biological_relevence_for_FM_feature/"

# 保存涉及相关性基因的所有连接
output_file1 <- paste0(output_dir, "relevant_protein_links_with_correlation.csv")
cat("正在保存涉及相关性基因的连接到：", output_file1, "\n")
fwrite(relevant_links_with_corr, output_file1, sep = ",")

# 保存两个基因都相关的连接
output_file2 <- paste0(output_dir, "both_genes_relevant_links_with_correlation.csv")
cat("正在保存两个基因都相关的连接到：", output_file2, "\n")
fwrite(both_genes_with_corr, output_file2, sep = ",")

# 创建网络分析（如果有足够的连接）
if(nrow(both_genes_with_corr) > 0) {
  cat("正在进行网络分析...\n")
  
  # 创建网络图
  edges <- both_genes_with_corr[, c(protein1_col, protein2_col), with = FALSE]
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 网络统计
  cat("网络统计：\n")
  cat("节点数量：", vcount(g), "\n")
  cat("边数量：", ecount(g), "\n")
  cat("网络密度：", round(edge_density(g), 4), "\n")
  cat("连通分量数量：", components(g)$no, "\n")
  
  # 计算节点度数
  node_degrees <- degree(g)
  cat("平均度数：", round(mean(node_degrees), 2), "\n")
  cat("最大度数：", max(node_degrees), "\n")
  
  # 找到度数最高的基因
  top_genes <- names(sort(node_degrees, decreasing = TRUE))[1:min(10, length(node_degrees))]
  cat("度数最高的基因（前10个）：\n")
  for(i in seq_along(top_genes)) {
    gene <- top_genes[i]
    corr <- gene_to_corr[gene]
    cat(sprintf("%d. %s (度数: %d, 相关性: %.4f)\n", 
                i, gene, node_degrees[gene], ifelse(is.na(corr), 0, corr)))
  }
  
  # 保存网络信息
  network_summary <- data.frame(
    Gene = names(node_degrees),
    Degree = node_degrees,
    Correlation = gene_to_corr[names(node_degrees)]
  )
  network_summary <- network_summary[order(network_summary$Degree, decreasing = TRUE), ]
  
  output_file3 <- paste0(output_dir, "network_gene_summary.csv")
  cat("正在保存网络基因摘要到：", output_file3, "\n")
  fwrite(network_summary, output_file3, sep = ",")
}

# 生成最终报告
cat("\n=== 分析完成报告 ===\n")
cat("1. 输入文件：\n")
cat("   - 相关性文件：", correlation_file, "\n")
cat("   - 连接文件：", converted_links_file, "\n")
cat("2. 基因统计：\n")
cat("   - 相关性基因数量：", length(correlation_genes), "\n")
cat("   - 连接文件基因数量：", length(protein_genes), "\n")
cat("   - 共同基因数量：", length(common_genes), "\n")
cat("3. 连接统计：\n")
cat("   - 原始连接数量：", nrow(protein_links), "\n")
cat("   - 涉及相关性基因的连接：", nrow(relevant_links), "\n")
cat("   - 两个基因都相关的连接：", nrow(both_genes_with_corr), "\n")
cat("4. 输出文件：\n")
cat("   - ", output_file1, "\n")
cat("   - ", output_file2, "\n")
if(exists("output_file3")) {
  cat("   - ", output_file3, "\n")
}

cat("\n脚本执行完成！\n")
