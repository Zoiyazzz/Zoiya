# 🧬 GSE32323｜基于DEG分析的生物知识图谱构建项目

本项目基于 GEO 公共数据集 GSE32323，通过差异表达分析与富集注释，构建了一个结构化的知识图谱 DEG → GO / KEGG，实现可视化、可查询的功能解读流程，适用于科研辅助、功能注释平台开发与知识挖掘应用。
本项目结合图数据库技术，将基因表达数据与功能富集注释转化为结构化知识图谱，实现更直观、更灵活的数据解读方式，适用于科研支持、功能注释平台开发与图挖掘应用场景。

---

## 🪧 项目概览

- **数据来源**：GSE32323（人类结直肠癌患者配对样本 + 5种癌细胞系）
- **测序平台**：Affymetrix GeneChip HG-U133 Plus 2.0（微阵列探针）
- **分析方式**：limma 差异分析 + clusterProfiler 富集注释
- **图谱平台**：Neo4j Aura（可视化 + Cypher 语义查询）
- **数据结构**：包含 Gene、GO_term、KEGG_pathway 三类节点，及其之间的注释边关系
- **应用目标**：
  - 功能富集关系结构化展示
  - 差异基因功能聚类与追踪
  - 支持标准化图谱构建流程与复用

---

## 📈 分析流程

```mermaid
graph LR
  A[表达矩阵下载] --> B[差异分析（limma）]
  B --> C[提取DEG（logFC + padj筛选）]
  C --> D[富集分析（GO / KEGG）]
  D --> E[构建图谱节点与边CSV]
  E --> F[图谱构建与交互查询]
```

## 🌐 图谱可视化（Neo4j Aura）
- **• DEG → GO term 局部网络图**
  <img width="890" height="581" alt="image" src="https://github.com/user-attachments/assets/ce5b38a1-82fe-4592-9df6-1c649ddf85b3" />
- **• DEG → KEGG path 局部网络图**
  <img width="892" height="587" alt="image" src="https://github.com/user-attachments/assets/6951aa29-2d84-4ed4-8463-2aa10e8ba253" />

- **🧠 查询示例与生物学解释**
- **1. 查询某功能对应 DEG**
```cypher
MATCH(g:Nodes human)-[:go_ human]->(go:Nodes human)
WHERE go.id ="G0:0015980"
RETURN g.id, g.logFC, go.description, g.`p.adjust!
0RDER BY abs(g.logFC) DESC  
```
<img width="892" height="625" alt="image" src="https://github.com/user-attachments/assets/4e5813a9-4da8-4c7e-8727-8e6df422a1ca" />

*结果示例：AGL, NDUFA1 等参与 “能量代谢” 的 DEG，可能反映样本代谢活动升高。*

- **2 查询某基因参与的GO功能**
```cypher
MATCH (g:Nodes_human {id: "ADHFE1"}) -I: go_human] ->(go:Nodes_human)
RETURN go.id, go description, go adjustedPValue
```
<img width="1790" height="1230" alt="image" src="https://github.com/user-attachments/assets/f339bde8-e551-459a-b948-df9e747c6161" />

*结果示例：ADHFE1在多个 分解代谢类 GO term 中出现，说明其核心功能与分解代谢通路密切相关。*

- **3 富集到 catabolic的所有 DEG 查询**
```cypher
MATCH (g:Nodes_human) - l:go_human, - (go: Nodes_human)
WHERE go description CONTAINS "catabolic"
RETURN g.id, g.logFC
```
<img width="1820" height="1182" alt="image" src="https://github.com/user-attachments/assets/d33a8421-2918-427c-b125-28ec248c9b31" />

*结果示例：可用于聚焦 catabolic相关通路，识别潜在机制关键因子。*

## 🎯项目价值总结

- ✅ 将富集分析从静态表格结构转化为图谱，可视化 DEG-功能关系
- ✅ 支持灵活语义查询：某功能→对应基因？某基因→哪些功能？
- ✅ 可复用于多个 bulk RNA-seq 项目，标准化知识抽取流程
- ✅ 适合作为知识图谱工程、图数据分析、生信挖掘等岗位展示项目

## 文件说明
- `DEG_results.csv`：原始差异表达结果
- `go_edges.csv`：基因与GO的注释关系
- `kegg_edges.csv`：基因与KEGG的注释关系
- `*_patch.csv`：用于补全节点属性的描述与p值
- `r_analysis.R`：R分析脚本
- `cypher_templates.md`：cypher 查询模版合集
