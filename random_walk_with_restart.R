# R脚本：Random Walk with Restart (RWR) 算法计算基因网络影响力
# 功能：基于PPI网络和初始相关性计算基因的网络影响力分数
# 算法：RWR with convergence detection using Euclidean distance
# 作者：AI助手
# 日期：2025-07-22

# 加载必要的包
library(data.table)
library(Matrix)
library(igraph)
library(dplyr)

# 设置文件路径
ppi_network_file <- "D:/biological_relevence_for_FM_feature/complete_network_with_self_loops.csv"
initial_correlation_file <- "D:/biological_relevence_for_FM_feature/tryfile_FTV_T0_corr_res_spear_train.csv"

# RWR算法参数设置
restart_probability <- 0.15  # 重启概率 (通常设置为0.1-0.2)
convergence_threshold <- 1e-6  # 收敛阈值
max_iterations <- 1000  # 最大迭代次数
distance_metric <- "euclidean"  # 距离度量方式: "euclidean" 或 "max_abs"

cat("=== Random Walk with Restart (RWR) 算法 ===\n")
cat("参数设置：\n")
cat("- 重启概率：", restart_probability, "\n")
cat("- 收敛阈值：", convergence_threshold, "\n")
cat("- 最大迭代次数：", max_iterations, "\n")
cat("- 距离度量：", distance_metric, "\n\n")

# 读取PPI网络文件
cat("正在读取PPI网络文件...\n")
ppi_network <- fread(ppi_network_file, header = TRUE)

cat("PPI网络文件结构：\n")
print(head(ppi_network))
cat("网络边数：", nrow(ppi_network), "\n")

# 获取网络中的所有基因
protein1_col <- names(ppi_network)[1]
protein2_col <- names(ppi_network)[2]
all_genes <- unique(c(ppi_network[[protein1_col]], ppi_network[[protein2_col]]))
all_genes <- all_genes[!is.na(all_genes)]  # 移除NA值
n_genes <- length(all_genes)

cat("网络中基因总数：", n_genes, "\n")

# 读取初始相关性文件
cat("正在读取初始相关性文件...\n")
initial_corr <- fread(initial_correlation_file, header = TRUE)

cat("相关性文件结构：\n")
print(head(initial_corr))
cat("相关性文件基因数：", nrow(initial_corr), "\n")

# 标准化列名
if(ncol(initial_corr) >= 2) {
  colnames(initial_corr)[1:2] <- c("Gene_Name", "Correlation")
}

# 构建转移概率矩阵
cat("正在构建转移概率矩阵...\n")

# 创建基因到索引的映射
gene_to_index <- setNames(1:n_genes, all_genes)
index_to_gene <- setNames(all_genes, 1:n_genes)

# 初始化邻接矩阵
adjacency_matrix <- Matrix(0, nrow = n_genes, ncol = n_genes, sparse = TRUE)
rownames(adjacency_matrix) <- all_genes
colnames(adjacency_matrix) <- all_genes

# 填充邻接矩阵
for(i in 1:nrow(ppi_network)) {
  gene1 <- ppi_network[[protein1_col]][i]
  gene2 <- ppi_network[[protein2_col]][i]
  
  if(!is.na(gene1) && !is.na(gene2) && gene1 %in% all_genes && gene2 %in% all_genes) {
    idx1 <- gene_to_index[gene1]
    idx2 <- gene_to_index[gene2]
    
    # 对于自连接，权重设为1；对于其他连接，也设为1（可以根据需要调整）
    adjacency_matrix[idx1, idx2] <- 1
    if(gene1 != gene2) {  # 避免重复设置自连接
      adjacency_matrix[idx2, idx1] <- 1
    }
  }
}

cat("邻接矩阵构建完成\n")
cat("非零元素数量：", nnzero(adjacency_matrix), "\n")

# 计算度数并构建转移概率矩阵
node_degrees <- rowSums(adjacency_matrix)
cat("度数统计：\n")
cat("- 最小度数：", min(node_degrees), "\n")
cat("- 最大度数：", max(node_degrees), "\n")
cat("- 平均度数：", round(mean(node_degrees), 2), "\n")

# 处理度数为0的节点（孤立节点）
isolated_nodes <- which(node_degrees == 0)
if(length(isolated_nodes) > 0) {
  cat("发现", length(isolated_nodes), "个孤立节点，为其添加自连接\n")
  for(idx in isolated_nodes) {
    adjacency_matrix[idx, idx] <- 1
  }
  node_degrees <- rowSums(adjacency_matrix)
}

# 构建转移概率矩阵 P
transition_matrix <- adjacency_matrix
for(i in 1:n_genes) {
  if(node_degrees[i] > 0) {
    transition_matrix[i, ] <- transition_matrix[i, ] / node_degrees[i]
  }
}

cat("转移概率矩阵构建完成\n")

# 构建初始概率向量
cat("正在构建初始概率向量...\n")

# 创建基因到相关性的映射
gene_to_corr <- setNames(initial_corr$Correlation, initial_corr$Gene_Name)

# 初始化概率向量
initial_prob <- numeric(n_genes)
names(initial_prob) <- all_genes

# 设置初始概率
for(i in 1:n_genes) {
  gene <- index_to_gene[i]
  if(gene %in% names(gene_to_corr)) {
    # 使用相关性的绝对值作为初始概率
    initial_prob[i] <- abs(gene_to_corr[gene])
  } else {
    # 对于没有相关性数据的基因，设置为很小的值
    initial_prob[i] <- 1e-10
  }
}

# 归一化初始概率向量
initial_prob <- initial_prob / sum(initial_prob)

cat("初始概率向量构建完成\n")
cat("有相关性数据的基因数：", sum(all_genes %in% names(gene_to_corr)), "\n")
cat("初始概率向量和：", sum(initial_prob), "\n")

# RWR算法实现
cat("\n=== 开始RWR算法迭代 ===\n")

# 初始化
current_prob <- initial_prob
iteration <- 0
converged <- FALSE
convergence_history <- numeric()

# 定义距离计算函数
calculate_distance <- function(vec1, vec2, metric = "euclidean") {
  if(metric == "euclidean") {
    return(sqrt(sum((vec1 - vec2)^2)))
  } else if(metric == "max_abs") {
    return(max(abs(vec1 - vec2)))
  } else {
    stop("不支持的距离度量方式")
  }
}

# 迭代过程
start_time <- Sys.time()

while(!converged && iteration < max_iterations) {
  iteration <- iteration + 1
  
  # 保存上一次的概率向量
  prev_prob <- current_prob
  
  # RWR更新公式: p(t+1) = (1-α) * P^T * p(t) + α * p(0)
  # 其中α是重启概率，P是转移矩阵，p(0)是初始概率向量
  current_prob <- (1 - restart_probability) * (t(transition_matrix) %*% current_prob) + 
                  restart_probability * initial_prob
  
  # 确保概率向量归一化
  current_prob <- as.numeric(current_prob) / sum(current_prob)
  
  # 计算收敛距离
  distance <- calculate_distance(current_prob, prev_prob, distance_metric)
  convergence_history <- c(convergence_history, distance)
  
  # 检查收敛
  if(distance < convergence_threshold) {
    converged <- TRUE
    cat("算法在第", iteration, "次迭代后收敛\n")
    cat("最终距离：", sprintf("%.2e", distance), "\n")
  }
  
  # 每100次迭代输出进度
  if(iteration %% 100 == 0) {
    cat("迭代", iteration, "次，当前距离：", sprintf("%.2e", distance), "\n")
  }
}

end_time <- Sys.time()
computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

if(!converged) {
  cat("警告：算法在", max_iterations, "次迭代后未收敛\n")
  cat("最终距离：", sprintf("%.2e", convergence_history[length(convergence_history)]), "\n")
}

cat("计算完成！\n")
cat("总迭代次数：", iteration, "\n")
cat("计算时间：", round(computation_time, 2), "秒\n")

# 整理结果
cat("\n=== 整理RWR结果 ===\n")

# 创建结果数据框
rwr_results <- data.frame(
  Gene = all_genes,
  RWR_Score = current_prob,
  Initial_Correlation = sapply(all_genes, function(g) {
    if(g %in% names(gene_to_corr)) gene_to_corr[g] else NA
  }),
  Node_Degree = node_degrees,
  stringsAsFactors = FALSE
)

# 按RWR分数排序
rwr_results <- rwr_results[order(rwr_results$RWR_Score, decreasing = TRUE), ]
rwr_results$Rank <- 1:nrow(rwr_results)

# 显示前20个基因
cat("RWR影响力分数前20的基因：\n")
print(head(rwr_results, 20))

# 统计信息
cat("\nRWR分数统计：\n")
cat("- 最小值：", sprintf("%.6e", min(rwr_results$RWR_Score)), "\n")
cat("- 最大值：", sprintf("%.6e", max(rwr_results$RWR_Score)), "\n")
cat("- 平均值：", sprintf("%.6e", mean(rwr_results$RWR_Score)), "\n")
cat("- 标准差：", sprintf("%.6e", sd(rwr_results$RWR_Score)), "\n")

# 保存结果
output_dir <- "D:/biological_relevence_for_FM_feature/"

# 保存RWR结果
rwr_output_file <- paste0(output_dir, "RWR_gene_influence_scores.csv")
cat("正在保存RWR结果到：", rwr_output_file, "\n")
fwrite(rwr_results, rwr_output_file, sep = ",")

# 保存收敛历史
convergence_output_file <- paste0(output_dir, "RWR_convergence_history.csv")
convergence_df <- data.frame(
  Iteration = 1:length(convergence_history),
  Distance = convergence_history
)
fwrite(convergence_df, convergence_output_file, sep = ",")

# 保存算法参数和统计信息
summary_file <- paste0(output_dir, "RWR_algorithm_summary.txt")
summary_info <- c(
  "=== Random Walk with Restart (RWR) 算法总结 ===",
  paste("执行时间：", Sys.time()),
  "",
  "算法参数：",
  paste("- 重启概率：", restart_probability),
  paste("- 收敛阈值：", convergence_threshold),
  paste("- 距离度量：", distance_metric),
  paste("- 最大迭代次数：", max_iterations),
  "",
  "网络统计：",
  paste("- 基因总数：", n_genes),
  paste("- 网络边数：", nrow(ppi_network)),
  paste("- 平均度数：", round(mean(node_degrees), 2)),
  "",
  "算法结果：",
  paste("- 实际迭代次数：", iteration),
  paste("- 是否收敛：", ifelse(converged, "是", "否")),
  paste("- 最终距离：", sprintf("%.2e", convergence_history[length(convergence_history)])),
  paste("- 计算时间：", round(computation_time, 2), "秒"),
  "",
  "输出文件：",
  paste("- RWR结果：", basename(rwr_output_file)),
  paste("- 收敛历史：", basename(convergence_output_file)),
  paste("- 算法总结：", basename(summary_file))
)

writeLines(summary_info, summary_file)

cat("\n=== RWR算法执行完成 ===\n")
cat("所有结果已保存到指定目录\n")
cat("主要输出文件：\n")
cat("1.", rwr_output_file, "\n")
cat("2.", convergence_output_file, "\n")
cat("3.", summary_file, "\n")
