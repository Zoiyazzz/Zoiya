# DEG → GO/KEGG 知识图谱构建项目（GSE32323）

本项目基于 GEO 数据集 GSE32323，完成从差异表达分析（DEG）到富集注释、知识图谱构建的全流程。图谱可在 Neo4j 中可视化、语义化查询。

## 项目亮点
- 差异表达分析（limma）
- GO/KEGG 富集分析（clusterProfiler）
- 图数据库建模（Neo4j Aura）
- 数据结构标准化（CSV + patch）
- Cypher 查询模版集合

## 文件说明
- `DEG_results.csv`：原始差异表达结果
- `go_edges.csv`：基因与GO的注释关系
- `kegg_edges.csv`：基因与KEGG的注释关系
- `*_patch.csv`：用于补全节点属性的描述与p值
- `r_analysis.R`：R分析脚本
