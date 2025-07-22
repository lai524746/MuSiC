# R脚本：将protein连接文件中的ENSP名字转换为symbol ID
# 作者：AI助手
# 日期：2025-07-22

# 加载必要的包
library(data.table)
library(dplyr)

# 设置文件路径
info_file <- "D:/biological_relevence_for_FM_feature/9606.protein.info.v12.0.txt"
links_file <- "D:/biological_relevence_for_FM_feature/9606.protein.links.v12.0.txt"

# 读取protein信息文件（ENSP到symbol ID的映射）
cat("正在读取protein信息文件...\n")
# 先尝试自动检测分隔符
protein_info <- fread(info_file, header = TRUE)

# 如果列数不对，尝试不同的分隔符
if(ncol(protein_info) == 1) {
  cat("尝试使用制表符分隔...\n")
  protein_info <- fread(info_file, header = TRUE, sep = "\t")
}
if(ncol(protein_info) == 1) {
  cat("尝试使用空格分隔...\n")
  protein_info <- fread(info_file, header = TRUE, sep = " ")
}

# 检查文件结构
cat("Protein信息文件结构：\n")
print(head(protein_info))
cat("信息文件维度：", dim(protein_info), "\n")
cat("信息文件列名：", names(protein_info), "\n")

# 假设第一列是ENSP，第二列是symbol ID
# 如果列名不同，请根据实际情况调整
colnames(protein_info)[1:2] <- c("ENSP", "Symbol_ID")

# 检查并清理protein_info数据
cat("清理protein信息数据...\n")
protein_info <- protein_info[!is.na(protein_info$ENSP) & !is.na(protein_info$Symbol_ID)]
protein_info <- protein_info[protein_info$ENSP != "" & protein_info$Symbol_ID != ""]

# 创建ENSP到Symbol ID的映射字典
ensp_to_symbol <- setNames(protein_info$Symbol_ID, protein_info$ENSP)
cat("映射字典大小：", length(ensp_to_symbol), "\n")

# 显示映射字典的示例
cat("映射字典示例：\n")
print(head(ensp_to_symbol, 10))

# 读取protein连接文件
cat("正在读取protein连接文件...\n")
# 先尝试空格分隔
protein_links <- fread(links_file, header = TRUE, sep = " ")

# 如果列数不对，尝试其他分隔符
if(ncol(protein_links) == 1) {
  cat("尝试使用制表符分隔...\n")
  protein_links <- fread(links_file, header = TRUE, sep = "\t")
}
if(ncol(protein_links) == 1) {
  cat("尝试自动检测分隔符...\n")
  protein_links <- fread(links_file, header = TRUE)
}

# 检查连接文件结构
cat("Protein连接文件结构：\n")
print(head(protein_links))
cat("连接文件维度：", dim(protein_links), "\n")
cat("连接文件列名：", names(protein_links), "\n")

# 备份原始数据
protein_links_original <- copy(protein_links)

# 假设连接文件的前两列是protein1和protein2的ENSP ID
# 根据实际文件结构调整列名
protein1_col <- names(protein_links)[1]
protein2_col <- names(protein_links)[2]

cat("要转换的列名：", protein1_col, "和", protein2_col, "\n")

# 检查连接文件中的ID格式
cat("连接文件中protein1的示例ID：\n")
print(head(unique(protein_links[[protein1_col]]), 10))
cat("连接文件中protein2的示例ID：\n")
print(head(unique(protein_links[[protein2_col]]), 10))

# 检查ID格式是否匹配
sample_ids_links <- head(unique(protein_links[[protein1_col]]), 5)
sample_ids_info <- head(unique(protein_info$ENSP), 5)
cat("连接文件ID示例：", sample_ids_links, "\n")
cat("信息文件ID示例：", sample_ids_info, "\n")

# 检查有多少ID能够匹配
protein1_ids <- unique(protein_links[[protein1_col]])
protein2_ids <- unique(protein_links[[protein2_col]])
all_link_ids <- unique(c(protein1_ids, protein2_ids))

match_count <- sum(all_link_ids %in% names(ensp_to_symbol))
cat("连接文件中唯一ID总数：", length(all_link_ids), "\n")
cat("能够匹配的ID数量：", match_count, "\n")
cat("匹配率：", round(match_count/length(all_link_ids)*100, 2), "%\n")

# 安全的转换函数
safe_convert <- function(ids, mapping_dict) {
  result <- character(length(ids))
  for(i in seq_along(ids)) {
    if(ids[i] %in% names(mapping_dict)) {
      result[i] <- mapping_dict[ids[i]]
    } else {
      result[i] <- NA_character_
    }
  }
  return(result)
}

cat("正在转换protein1列的ENSP ID到Symbol ID...\n")
# 使用安全转换函数
protein_links[[protein1_col]] <- safe_convert(protein_links[[protein1_col]], ensp_to_symbol)

cat("正在转换protein2列的ENSP ID到Symbol ID...\n")
# 使用安全转换函数
protein_links[[protein2_col]] <- safe_convert(protein_links[[protein2_col]], ensp_to_symbol)

# 检查转换结果
cat("转换后的连接文件结构：\n")
print(head(protein_links))

# 检查是否有未能转换的ENSP ID（显示为NA）
na_count_protein1 <- sum(is.na(protein_links[[protein1_col]]))
na_count_protein2 <- sum(is.na(protein_links[[protein2_col]]))

cat("转换统计：\n")
cat("Protein1列中未能转换的ENSP ID数量：", na_count_protein1, "\n")
cat("Protein2列中未能转换的ENSP ID数量：", na_count_protein2, "\n")
cat("总行数：", nrow(protein_links), "\n")
cat("Protein1转换成功率：", round((1 - na_count_protein1/nrow(protein_links))*100, 2), "%\n")
cat("Protein2转换成功率：", round((1 - na_count_protein2/nrow(protein_links))*100, 2), "%\n")

# 显示一些未能转换的ID示例（如果有的话）
if(na_count_protein1 > 0) {
  failed_ids_1 <- protein_links_original[[protein1_col]][is.na(protein_links[[protein1_col]])]
  cat("未能转换的Protein1 ID示例：\n")
  print(head(unique(failed_ids_1), 10))
}

if(na_count_protein2 > 0) {
  failed_ids_2 <- protein_links_original[[protein2_col]][is.na(protein_links[[protein2_col]])]
  cat("未能转换的Protein2 ID示例：\n")
  print(head(unique(failed_ids_2), 10))
}

# 如果需要，可以移除包含NA的行
if(na_count_protein1 > 0 || na_count_protein2 > 0) {
  cat("是否要移除包含未转换ID的行？(y/n)\n")
  # 在实际使用中，您可以取消注释下面的行来交互式选择
  # user_choice <- readline()
  # 这里我们默认保留所有行，包括NA
  user_choice <- "n"
  
  if(tolower(user_choice) == "y") {
    protein_links_clean <- protein_links[!is.na(protein_links[[protein1_col]]) & 
                                       !is.na(protein_links[[protein2_col]])]
    cat("清理后的文件维度：", dim(protein_links_clean), "\n")
  } else {
    protein_links_clean <- protein_links
  }
} else {
  protein_links_clean <- protein_links
  cat("所有ENSP ID都成功转换为Symbol ID！\n")
}

# 保存转换后的文件
output_file <- "D:/biological_relevence_for_FM_feature/9606.protein.links.v12.0_symbol_converted.txt"
cat("正在保存转换后的文件到：", output_file, "\n")

fwrite(protein_links_clean, output_file, sep = "\t", quote = FALSE)

cat("转换完成！\n")
cat("原始连接文件行数：", nrow(protein_links_original), "\n")
cat("转换后文件行数：", nrow(protein_links_clean), "\n")

# 显示转换前后的对比示例
cat("\n转换前后对比示例：\n")
comparison_sample <- data.frame(
  原始_Protein1 = head(protein_links_original[[protein1_col]], 5),
  转换后_Protein1 = head(protein_links_clean[[protein1_col]], 5),
  原始_Protein2 = head(protein_links_original[[protein2_col]], 5),
  转换后_Protein2 = head(protein_links_clean[[protein2_col]], 5)
)
print(comparison_sample)

# 可选：创建转换统计报告
cat("\n生成转换统计报告...\n")
conversion_stats <- list(
  原始文件行数 = nrow(protein_links_original),
  转换后文件行数 = nrow(protein_links_clean),
  映射字典大小 = length(ensp_to_symbol),
  未转换的Protein1数量 = na_count_protein1,
  未转换的Protein2数量 = na_count_protein2,
  转换成功率_Protein1 = round((1 - na_count_protein1/nrow(protein_links_original)) * 100, 2),
  转换成功率_Protein2 = round((1 - na_count_protein2/nrow(protein_links_original)) * 100, 2)
)

print(conversion_stats)

cat("\n脚本执行完成！\n")
