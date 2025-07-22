# R脚本：基于基因表达数据构建蛋白质相互作用网络
# 功能：从基因表达文件中提取基因名，并在蛋白质相互作用网络中找到相应连接
# 特色：防止基因名被Excel误读为日期格式（如SEPT4等）
# 作者：AI助手
# 日期：2025-07-22

# 加载必要的包
library(data.table)
library(dplyr)
library(igraph)

# 设置文件路径
gene_expression_file <- "D:/biological_relevence_for_FM_feature/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt"
converted_links_file <- "D:/biological_relevence_for_FM_feature/9606.protein.links.v12.0_symbol_converted.txt"

# 读取基因表达文件（防止基因名被转换为日期）
cat("正在读取基因表达文件...\n")
cat("文件路径：", gene_expression_file, "\n")

# 使用多种方法防止基因名被转换为日期
# 方法1：使用fread并指定colClasses
gene_expression <- tryCatch({
  fread(gene_expression_file, header = TRUE, stringsAsFactors = FALSE,
        colClasses = list(character = 1))  # 第一列强制为字符型
}, error = function(e) {
  cat("fread失败，尝试read.table...\n")
  # 方法2：使用read.table并禁用字符串转换
  read.table(gene_expression_file, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE,
             row.names = 1, colClasses = c("character", rep("numeric", 100)))
})

# 如果上述方法都失败，尝试手动读取
if(is.null(gene_expression) || nrow(gene_expression) == 0) {
  cat("尝试手动读取文件...\n")
  # 先读取第一行获取列数
  first_line <- readLines(gene_expression_file, n = 1)
  ncols <- length(strsplit(first_line, "\t")[[1]])

  # 设置所有列的类型，第一列为字符型
  col_classes <- c("character", rep("numeric", ncols - 1))

  gene_expression <- read.table(gene_expression_file, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               colClasses = col_classes)
}

# 检查文件结构
cat("基因表达文件结构：\n")
print(head(gene_expression[, 1:min(5, ncol(gene_expression))]))  # 只显示前5列
cat("基因表达文件维度：", dim(gene_expression), "\n")
cat("基因表达文件列名（前10个）：", head(names(gene_expression), 10), "\n")

# 提取基因名（行名或第一列）
if(is.null(rownames(gene_expression)) || all(rownames(gene_expression) == as.character(1:nrow(gene_expression)))) {
  # 如果没有行名或行名是默认数字，使用第一列作为基因名
  gene_names <- gene_expression[[1]]
  cat("使用第一列作为基因名\n")
} else {
  # 使用行名作为基因名
  gene_names <- rownames(gene_expression)
  cat("使用行名作为基因名\n")
}

# 清理基因名：移除空值、空字符串和重复项
gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]
gene_names <- unique(gene_names)

cat("提取的基因数量：", length(gene_names), "\n")
cat("基因名示例：\n")
print(head(gene_names, 15))

# 检查是否有可能被误读为日期的基因名
date_like_genes <- gene_names[grepl("^[0-9]{1,2}[-/][0-9]{1,2}$|^[A-Z]{3}[0-9]{1,2}$|^[0-9]{1,2}-[A-Z]{3}$", gene_names)]
if(length(date_like_genes) > 0) {
  cat("发现可能被误读为日期的基因名：\n")
  print(head(date_like_genes, 10))
  cat("总数：", length(date_like_genes), "\n")
} else {
  cat("未发现明显的日期格式基因名，读取成功！\n")
}

# 创建基因列表数据框（模拟相关性数据结构以保持代码兼容性）
gene_list_df <- data.frame(
  Gene_Name = gene_names,
  Present_in_Expression = 1,  # 标记基因在表达文件中存在
  stringsAsFactors = FALSE
)

cat("基因列表数据框维度：", dim(gene_list_df), "\n")

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

# 获取基因表达文件中的基因列表
expression_genes <- unique(gene_list_df$Gene_Name)
cat("基因表达文件中的基因数量：", length(expression_genes), "\n")

# 显示基因示例
cat("基因名示例：\n")
print(head(expression_genes, 10))

# 获取protein连接文件中的所有基因
protein_genes <- unique(c(protein_links[[protein1_col]], protein_links[[protein2_col]]))
# 移除NA值
protein_genes <- protein_genes[!is.na(protein_genes)]
cat("Protein连接文件中的基因数量：", length(protein_genes), "\n")

# 找到两个文件中共同的基因
common_genes <- intersect(expression_genes, protein_genes)
cat("共同基因数量：", length(common_genes), "\n")
cat("覆盖率（表达基因在连接文件中的比例）：",
    round(length(common_genes)/length(expression_genes)*100, 2), "%\n")

# 显示共同基因示例
if(length(common_genes) > 0) {
  cat("共同基因示例：\n")
  print(head(common_genes, 10))
} else {
  cat("警告：没有找到共同的基因！请检查基因名格式是否一致。\n")
  cat("表达文件基因示例：", head(expression_genes, 5), "\n")
  cat("连接文件基因示例：", head(protein_genes, 5), "\n")
}

# 筛选出涉及表达基因的protein连接
cat("正在筛选涉及表达基因的protein连接...\n")
relevant_links <- protein_links[
  (protein_links[[protein1_col]] %in% expression_genes) |
  (protein_links[[protein2_col]] %in% expression_genes)
]

cat("筛选后的连接数量：", nrow(relevant_links), "\n")
cat("原始连接数量：", nrow(protein_links), "\n")
cat("筛选比例：", round(nrow(relevant_links)/nrow(protein_links)*100, 2), "%\n")

# 进一步筛选：只保留两个基因都在表达列表中的连接
cat("正在筛选两个基因都在表达列表中的连接...\n")
both_genes_relevant <- relevant_links[
  (relevant_links[[protein1_col]] %in% expression_genes) &
  (relevant_links[[protein2_col]] %in% expression_genes)
]

cat("两个基因都在表达文件中的连接数量：", nrow(both_genes_relevant), "\n")

# 为连接添加表达信息
cat("正在为连接添加基因表达信息...\n")

# 创建基因到表达状态的映射（所有基因都标记为存在）
gene_to_expression <- setNames(rep(1, length(expression_genes)), expression_genes)

# 为relevant_links添加表达信息
relevant_links_with_info <- relevant_links
relevant_links_with_info$Gene1_In_Expression <- gene_to_expression[relevant_links_with_info[[protein1_col]]]
relevant_links_with_info$Gene2_In_Expression <- gene_to_expression[relevant_links_with_info[[protein2_col]]]

# 为both_genes_relevant添加表达信息
both_genes_with_info <- both_genes_relevant
both_genes_with_info$Gene1_In_Expression <- gene_to_expression[both_genes_with_info[[protein1_col]]]
both_genes_with_info$Gene2_In_Expression <- gene_to_expression[both_genes_with_info[[protein2_col]]]

# 计算表达统计
if(nrow(both_genes_with_info) > 0) {
  cat("表达统计（两个基因都在表达文件中的连接）：\n")
  cat("连接中基因1在表达文件中的数量：", sum(!is.na(both_genes_with_info$Gene1_In_Expression)), "\n")
  cat("连接中基因2在表达文件中的数量：", sum(!is.na(both_genes_with_info$Gene2_In_Expression)), "\n")
  cat("两个基因都在表达文件中的连接比例：100%\n")
}

# 保存结果文件
output_dir <- "D:/biological_relevence_for_FM_feature/"

# 保存涉及表达基因的所有连接
output_file1 <- paste0(output_dir, "relevant_protein_links_with_expression_genes.csv")
cat("正在保存涉及表达基因的连接到：", output_file1, "\n")
fwrite(relevant_links_with_info, output_file1, sep = ",")

# 保存两个基因都在表达文件中的连接
output_file2 <- paste0(output_dir, "both_genes_in_expression_links.csv")
cat("正在保存两个基因都在表达文件中的连接到：", output_file2, "\n")
fwrite(both_genes_with_info, output_file2, sep = ",")

# 创建网络分析（如果有足够的连接）
if(nrow(both_genes_with_info) > 0) {
  cat("正在进行网络分析...\n")

  # 创建网络图
  edges <- both_genes_with_info[, c(protein1_col, protein2_col), with = FALSE]
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
    in_expression <- gene_to_expression[gene]
    cat(sprintf("%d. %s (度数: %d, 在表达文件中: %s)\n",
                i, gene, node_degrees[gene], ifelse(is.na(in_expression), "否", "是")))
  }

  # 保存网络信息
  network_summary <- data.frame(
    Gene = names(node_degrees),
    Degree = node_degrees,
    In_Expression_File = gene_to_expression[names(node_degrees)],
    stringsAsFactors = FALSE
  )
  network_summary <- network_summary[order(network_summary$Degree, decreasing = TRUE), ]

  output_file3 <- paste0(output_dir, "network_gene_summary_expression.csv")
  cat("正在保存网络基因摘要到：", output_file3, "\n")
  fwrite(network_summary, output_file3, sep = ",")
}

# 生成最终报告
cat("\n=== 分析完成报告 ===\n")
cat("1. 输入文件：\n")
cat("   - 基因表达文件：", gene_expression_file, "\n")
cat("   - 连接文件：", converted_links_file, "\n")
cat("2. 基因统计：\n")
cat("   - 表达文件基因数量：", length(expression_genes), "\n")
cat("   - 连接文件基因数量：", length(protein_genes), "\n")
cat("   - 共同基因数量：", length(common_genes), "\n")
cat("   - 覆盖率：", round(length(common_genes)/length(expression_genes)*100, 2), "%\n")
cat("3. 连接统计：\n")
cat("   - 原始连接数量：", nrow(protein_links), "\n")
cat("   - 涉及表达基因的连接：", nrow(relevant_links), "\n")
cat("   - 两个基因都在表达文件中的连接：", nrow(both_genes_with_info), "\n")
cat("4. 输出文件：\n")
cat("   - ", output_file1, "\n")
cat("   - ", output_file2, "\n")
if(exists("output_file3")) {
  cat("   - ", output_file3, "\n")
}

# 保存可能被误读为日期的基因名列表
if(length(date_like_genes) > 0) {
  date_genes_file <- paste0(output_dir, "potential_date_format_genes.txt")
  writeLines(date_like_genes, date_genes_file)
  cat("   - ", date_genes_file, "\n")
}

cat("5. 基因名保护：\n")
if(length(date_like_genes) > 0) {
  cat("   - 检测到", length(date_like_genes), "个可能被误读为日期的基因名\n")
  cat("   - 已使用字符串保护措施防止自动转换\n")
  cat("   - 详细列表已保存到 potential_date_format_genes.txt\n")
} else {
  cat("   - 未检测到明显的日期格式基因名\n")
}

cat("\n=== 使用建议 ===\n")
cat("1. 在Excel中打开CSV文件时，建议：\n")
cat("   - 使用'数据' -> '从文本/CSV导入'功能\n")
cat("   - 将基因名列设置为'文本'格式\n")
cat("   - 避免直接双击打开CSV文件\n")
cat("2. 如果发现基因名仍被转换为日期，请检查 potential_date_format_genes.txt 文件\n")
cat("3. 建议在后续分析中继续使用R或Python等编程工具处理数据\n")

cat("\n脚本执行完成！\n")
