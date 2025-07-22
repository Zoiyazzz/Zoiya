# 📖 Cypher 查询语句模板合集（DEG-GO-KEGG 图谱）

## **📌 节点基础查询**

**查询所有节点类型（nodeType）**

```sql
#Cypher
MATCH (n)
RETURN DISTINCT n.nodeType
```

**查询所有关系类型**

```sql
#Cypher
MATCH ()-[r]->()
RETURN DISTINCT type(r)
```

# **🧬 DEG 查询相关**

**查询所有上调 DEG（logFC > 0）**

```sql
#Cypher
MATCH (g:Nodes_human)
WHERE g.log2FoldChange > 0
RETURN g.id, g.log2FoldChange, g.pvalue
ORDER BY g.log2FoldChange DESC
```

**查询所有下调 DEG（logFC < 0）**

```sql
#Cypher
MATCH (g:Nodes_human)
WHERE g.log2FoldChange < 0
RETURN g.id, g.log2FoldChange, g.pvalue
ORDER BY g.log2FoldChange ASC
```

**查询 DEG 的 logFC 最大 / 最小值**

```sql
#Cypher
MATCH (g:Nodes_human)
RETURN max(g.log2FoldChange) AS max_logFC, min(g.log2FoldChange) AS min_logFC
```

# **🔍 GO term 功能查询**

**某 GO term 富集到的 DEG**

```sql
#Cypher
MATCH (g:Nodes_human)-[:go_human]->(go:Nodes_human)
WHERE go.id = "GO:0044282"
RETURN g.id, g.log2FoldChange, g.pvalue
```

**查询包含关键词的 GO term（模糊搜索）**

```sql
#Cypher
MATCH (go:Nodes_human)
WHERE go.nodeType = "GO_term" AND go.description CONTAINS "catabolic"
RETURN [go.id](http://go.id/), go.description, go.adjustedPValue
```

**查询某基因参与的所有 GO term**

```sql
#Cypher
MATCH (g:Nodes_human {id: "ADHFE1"})-[:go_human]->(go:Nodes_human)
RETURN go.id, go.description, go.adjustedPValue
```

# **🧠 KEGG 通路查询**

**某通路富集的 DEG**

```sql
#Cypher
MATCH (g:Nodes_human)-[:kegg_human]->(k:Nodes_human)
WHERE k.id = "hsa04110"
RETURN g.id, g.log2FoldChange, g.pvalue
```

**查询所有 KEGG 通路及其 DEG 数量**

```sql
#Cypher
MATCH (g:Nodes_human)-[:kegg_human]->(k:Nodes_human)
RETURN k.id, k.description, count(g) AS degCount
ORDER BY degCount DESC
```

# **📊 综合查询：GO + KEGG 交叉分析**

**同时查询 DEG 对应的 GO 和 KEGG 注释**

```sql
#Cypher
MATCH (g:Nodes_human)-[:go_human]->(go:Nodes_human),
      (g)-[:kegg_human]->(k:Nodes_human)
RETURN g.id, go.description AS GO_desc, k.description AS KEGG_desc
```

# **🌐 图谱可视化结构用语句**

**返回某 GO term 与 DEG 的放射结构**

```sql
#Cypher
MATCH (g:Nodes_human)-[r:go_human]->(go:Nodes_human)
WHERE go.id = "GO:0015980"
RETURN g, r, go
```

**某 KEGG 通路的 DEG 放射结构图**

```sql
#Cypher
MATCH (g:Nodes_human)-[r:kegg_human]->(k:Nodes_human)
WHERE k.id = "hsa04110"
RETURN g, r, k
```

# **🎯 统计类查询**

**上调 / 下调 DEG 数量**

```sql
#Cypher
MATCH (g:Nodes_human)
WHERE exists(g.log2FoldChange)
RETURN
  count(CASE WHEN g.log2FoldChange > 0 THEN 1 END) AS upregulated,
  count(CASE WHEN g.log2FoldChange < 0 THEN 1 END) AS downregulated
```

**某 GO term 中上下调基因数量**

```sql
#Cypher
MATCH (g:Nodes_human)-[:go_human]->(go:Nodes_human)
WHERE go.id = "GO:0015980"
RETURN
  count(CASE WHEN g.log2FoldChange > 0 THEN 1 END) AS up,
  count(CASE WHEN g.log2FoldChange < 0 THEN 1 END) AS down
```