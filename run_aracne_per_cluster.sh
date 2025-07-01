#!/bin/bash

# ======== 配置部分 ========
TF_LIST=../tfs.txt
PVAL=1E-8
MEM=5G
ARACNE_JAR=aracne.jar

# 要处理的输入表达矩阵文件名（不带.tsv 后缀）
FILES=("venetoclax-r1-clusts_clust-1-metaCells" "venetoclax-r1-clusts_clust-2-metaCells")

# ======== 主循环 ========
for BASE in "${FILES[@]}"
do
  INPUT="../${BASE}.tsv"
  OUTPUT_DIR="../${BASE}_aracne"

  echo ">>> Processing $INPUT"

  # 创建输出目录
  mkdir -p "$OUTPUT_DIR"

  # Step 1: calculateThreshold
  java -Xmx$MEM -jar $ARACNE_JAR \
    -e "$INPUT" \
    -o "$OUTPUT_DIR" \
    --tfs "$TF_LIST" \
    --pvalue $PVAL \
    --seed 1 \
    --calculateThreshold

  # Step 2: 100 bootstraps
  for i in {1..100}
  do
    echo "  Bootstrap $i for $BASE"
    java -Xmx$MEM -jar $ARACNE_JAR \
      -e "$INPUT" \
      -o "$OUTPUT_DIR" \
      --tfs "$TF_LIST" \
      --pvalue $PVAL \
      --seed $i
  done

  # Step 3: consolidate
  java -Xmx$MEM -jar $ARACNE_JAR \
    -o "$OUTPUT_DIR" \
    --consolidate

  echo ">>> Done with $BASE"
done
