# GFF / GTF 格式

**GFF**（General Feature Format）與 **GTF**（Gene Transfer Format）是兩種描述基因組特徵的格式，均由 **9 個欄位**組成。前 8 欄完全相同，差異在於第 9 欄的寫法。

> GTF 本質上是 GFF 第 2 版（GFF2）的衍生格式。

---

### 欄位說明

| # | 欄位名稱 | 說明 |
|---|----------|------|
| 1 | `seqname` | 序列名稱，通常是染色體或 scaffold。染色體名稱可帶或不帶 `chr` 前綴，必須與 Ensembl 中使用的名稱一致 |
| 2 | `source` | 資料來源或 annotation pipeline 名稱，例如 `HAVANA`、`Ensembl` |
| 3 | `feature` | 特徵類型，例如 `gene`、`transcript`、`exon`、`CDS` |
| 4 | `start` | 特徵的起始位置（1-based） |
| 5 | `end` | 特徵的結束位置（1-based，包含） |
| 6 | `score` | 可信度分數或統計分數；若無則填 `.` |
| 7 | `strand` | 所在 DNA 鏈：`+` 為正鏈，`-` 為負鏈 |
| 8 | `frame` | 僅用於 CDS，表示 reading frame 從第幾個 base 開始讀 codon（`0`、`1`、`2`） |
| 9 | `attribute` | 附加資訊欄，詳見下方說明 |

---

### 第 9 欄：attribute（附加資訊）

這是最複雜也最重要的一欄。GTF 與 GFF3 在此欄的**格式**與**設計邏輯**都不同。

#### 格式差異

| | GTF（GFF2） | GFF3 |
|-|-------------|------|
| 分隔符 | `key "value";` 空白分隔 key 與 value | `key=value;` 等號分隔 |
| 範例 | `gene_id "ENSG00000160072"; gene_name "ATAD3B";` | `ID=gene:ENSG00000160072;Name=ATAD3B;` |

#### 階層關係的表達方式

描述 gene → transcript → exon 的方式不同：

- **GTF**：每一行都重複寫 `gene_id` 和 `transcript_id`，靠這兩個欄位隱式連結上下層
- **GFF3**：用 `ID` 標記自己，用 `Parent` 明確指向上層 feature，結構像樹狀圖

```
# GFF3 範例
gene    → ID=gene:ENSG00000160072;Name=ATAD3B
mRNA    → ID=transcript:ENST00000356607;Parent=gene:ENSG00000160072
exon    → Parent=transcript:ENST00000356607
```

---

#### GTF 常見 attribute 欄位

| attribute | 說明 |
|-----------|------|
| `gene_id` | Ensembl gene ID，例如 `ENSG00000160072` |
| `gene_name` | 常用基因名稱，例如 `ATAD3B` |
| `gene_version` | gene ID 的版本號 |
| `gene_biotype` | gene 類型，例如 `protein_coding`、`lncRNA` |
| `transcript_id` | Ensembl transcript ID，例如 `ENST00000...` |
| `transcript_name` | transcript 名稱 |
| `transcript_biotype` | transcript 類型 |
| `exon_number` | 第幾個 exon |
| `exon_id` | Ensembl exon ID |

#### GFF3 標準 attribute 欄位

| attribute | 說明 |
|-----------|------|
| `ID` | 該 feature 的唯一識別碼，在整個檔案中必須唯一 |
| `Name` | 顯示用名稱，例如基因名稱 `ATAD3B` |
| `Parent` | 指向上層 feature 的 `ID`，建立階層關係 |
| `Dbxref` | 外部資料庫連結，例如 `Dbxref=GeneID:23300` |
| `Ontology_term` | GO term 等功能分類 |

> **補充**：Ensembl 提供的 GFF3 檔案實際上也會附上 `gene_id`、`gene_name` 等欄位，與 GTF 共用部分 attribute 名稱，方便工具互通。

---

#### 搜尋特定基因（以 ATAD3B 為例）

GTF：
```bash
grep 'gene_name "ATAD3B"' annotation.gtf
grep 'gene_id "ENSG00000160072"' annotation.gtf
```

GFF3：
```bash
grep 'Name=ATAD3B' annotation.gff3
grep 'ID=gene:ENSG00000160072' annotation.gff3
```

因此如果純粹找基因完整的位子時，feature 為gene，attribute 為 gene_id 或 gene_name 就好：
1	ensembl_havana	gene	1471765	1497848	.	+	.	gene_id "ENSG00000160072"; gene_version "20"; gene_name "ATAD3B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";