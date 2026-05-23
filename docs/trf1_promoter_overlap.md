# TRF1 ChIP-seq promoter overlap analysis

這份流程用 TRF1 ChIP-seq peak 與 DEG gene list，檢查目標基因 TSS 上游指定距離內是否有 TRF1 binding。

核心邏輯固定是：

```text
gene annotation -> TSS -> upstream promoter window -> TRF1 peak overlap
```

沒有把 peak 直接指定到最近基因。

## 檔案

主要輸入與程式碼位置：

```text
data/GSM638201_TRF1-p0.05.bed.gz       # TRF1 ChIP-seq peak (SISSRs)
data/siTRF1_DEG_list_gene_name.txt     # DEG list (gene_id/gene_name/regulation)
data/gene_annotation.tsv               # GTF 轉成的 BED-style 註解
src/trf1_overlap/                      # pipeline 套件
scripts/trf1_promoter_overlap.py       # CLI wrapper
```

如果只有 `Homo_sapiens.GRCh38.108.gtf` 而沒有 `gene_annotation.tsv`，CLI 第一次跑會自動產生並快取：

```bash
conda run -n Bioinfomatics trf1-promoter-overlap \
  --peak data/GSM638201_TRF1-p0.05.bed.gz \
  --deg data/siTRF1_DEG_list_gene_name.txt \
  --annotation data/gene_annotation.tsv \
  --gtf Homo_sapiens.GRCh38.108.gtf \
  --outdir results
```

之後 `data/gene_annotation.tsv` 已存在，就可以不用 `--gtf`。

## Conda 環境

使用 `Bioinfomatics` conda environment。建立方式：

```bash
conda env create -f environment.yml
conda activate Bioinfomatics
```

`environment.yml` 已經包含 `pip install -e .`，所以 `trf1-promoter-overlap` CLI、`python -m trf1_overlap` 都可直接用。

主要依賴：

```text
python 3.12
pandas
matplotlib
pybedtools
bedtools
```

`pybedtools` 目前沒有 Python 3.14 conda build，所以環境鎖 Python 3.12。

## 正式執行

```bash
conda run -n Bioinfomatics trf1-promoter-overlap \
  --peak data/GSM638201_TRF1-p0.05.bed.gz \
  --deg data/siTRF1_DEG_list_gene_name.txt \
  --annotation data/gene_annotation.tsv \
  --outdir results
```

或 `make trf1`。

輸出：

```text
results/TRF1_cleaned_peaks.bed
results/TRF1_promoter_overlap_1000.tsv
results/TRF1_promoter_overlap_2000.tsv
results/TRF1_promoter_overlap_3000.tsv
results/TRF1_promoter_overlap_summary.tsv
results/TRF1_matched_genes_barplot.png
results/TRF1_matched_genes_stacked_barplot.png
results/report.txt
```

## 輸入格式

### Peak file

原始 peak 檔含 SISSRs header 與下列欄位：

```text
Chr cStart cEnd NumTags Fold p-value
```

腳本會清理成 BED-style peak：

```text
chr start end
```

同時保留 `Fold` 與 `p-value`，用於最後 overlap TSV 的 `peak_fold` 與 `peak_pvalue` 欄位。

### DEG list

需要欄位：

```text
gene_id
gene_name
regulation
```

`regulation` 只接受：

```text
UP
DOWN
```

### Gene annotation

`gene_annotation.tsv` 需要欄位：

```text
chr
gene_start
gene_end
gene_name
gene_id
strand
```

如果從 GTF 產生，腳本會把 GTF 的 1-based closed start 轉成 BED-style 0-based start，方便與 peak 用 bedtools intersect。

染色體名稱會正規化：

```text
1 -> chr1
MT -> chrM
M -> chrM
```

## Upstream promoter window 定義

TSS 規則：

```text
strand +: TSS = gene_start
strand -: TSS = gene_end
```

這裡的 1000bp、2000bp、3000bp 不是三個獨立 promoter，而是從同一個 TSS 往上游延伸的三個巢狀搜尋範圍。1000bp window 會被包含在 2000bp window 裡，2000bp window 會被包含在 3000bp window 裡。

Window 規則：

```text
strand +:
  1000bp = [TSS-1000, TSS]
  2000bp = [TSS-2000, TSS]
  3000bp = [TSS-3000, TSS]

strand -:
  1000bp = [TSS, TSS+1000]
  2000bp = [TSS, TSS+2000]
  3000bp = [TSS, TSS+3000]
```

window start 小於 0 時會修正為 0。輸出欄位仍命名為 `promoter_start` 與 `promoter_end`，代表這個 TSS upstream promoter window 的座標。

## 目前分析結果

正式分析結果：

```text
1000bp: total genes=274, matched genes=1, UP matched=1, DOWN matched=0
2000bp: total genes=274, matched genes=1, UP matched=1, DOWN matched=0
3000bp: total genes=274, matched genes=1, UP matched=1, DOWN matched=0
```

三個 upstream window threshold 的 matched gene 都是：

```text
MT-ND1    ENSG00000198888    UP
```

準確的結論寫法應該是：

```text
MT-ND1 的 TSS 上游 1000bp 內有 TRF1 peak；
因此在 TSS 上游 2000bp 與 3000bp window 中也會被包含。
```

不要寫成「TRF1 前 1000-3000bp 中有 MT-ND1」，因為 TRF1 在這裡是 ChIP-seq binding peak，不是一個用來定義上游區域的 gene。這個分析是 gene-centric：先定義每個 gene 的 TSS upstream window，再問 TRF1 peak 是否 overlap 這個 window。

這個結果是資料本身造成的，不是 report 彙整錯誤：

- `MT-ND1` 確實在 DEG list 裡，regulation 是 `UP`
- GTF 裡 `MT-ND1` 位於 mitochondrial chromosome，TSS 約在 `chrM:3306`
- TRF1 peak 檔有多個 `chrM` peaks
- 這些 `chrM` peaks 落在 `MT-ND1` 的 1kb、2kb、3kb upstream promoter window 內

例如 1000bp upstream window 有 3 個 peak overlap：

```text
MT-ND1 chrM upstream window 2306-3306 overlaps peaks:
chrM 2051-2471
chrM 2611-2731
chrM 2991-3151
```

因為 2000bp 與 3000bp upstream window 更大，會包含 1000bp window，所以同一個 gene 仍然會保留在結果中，只是額外納入更遠上游的 peaks。這不能解讀成同一個基因跨了三個不同地段，而是同一個 TSS 的搜尋距離被逐步放寬。

## Nuclear-only sanity check

如果分析目標只想看 nuclear genome promoter，而不想讓 mitochondrial peaks 進入結果，可以排除 `chrM`：

```bash
conda run -n Bioinfomatics trf1-promoter-overlap \
  --peak data/GSM638201_TRF1-p0.05.bed.gz \
  --deg data/siTRF1_DEG_list_gene_name.txt \
  --annotation data/gene_annotation.tsv \
  --outdir results_no_chrM \
  --exclude-chroms chrM
```

或 `make trf1-no-chrM`。

目前 sanity check 結果：

```text
1000bp: matched genes=0
2000bp: matched genes=0
3000bp: matched genes=0
```

也就是說，若排除 mitochondrial chromosome，這批 DEG genes 在 TSS 上游 1kb、2kb、3kb 範圍內沒有其他 TRF1 peak overlap。

## 常用參數

```text
--peak            TRF1 peak 檔，支援 .bed 或 .bed.gz
--deg             DEG gene list
--annotation      gene annotation TSV
--gtf             annotation TSV 不存在時，用此 GTF 產生 annotation
--outdir          輸出資料夾
--exclude-chroms  可選，逗號分隔要排除的染色體，例如 chrM 或 chrM,chrY
```

## 解讀注意事項

`MT-ND1` 是 mitochondrial gene。TRF1 是 telomere-associated factor，若研究問題聚焦在 nuclear promoter binding，建議把 `chrM` 結果視為 QC/sanity check 項目，並用 `--exclude-chroms chrM` 產生 nuclear-only 結果。

預設流程不會自動排除 `chrM`，因為這會改變原始輸入資料的分析範圍；排除與否應該由分析問題決定。
