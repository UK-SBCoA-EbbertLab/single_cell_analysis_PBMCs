## ExpressionMatrixProcessor

This Java processing script filters and transposes tab-delimited expression matrices (genes/isoforms × cells) to reduce file size for downstream supervised filtering.

It performs the following functions:
1. Removes low-depth cells  
2. Removes rarely expressed genes and isoforms  
3. Rounds numeric values to two decimal places (removes decimals if integer)  
4. Writes a filtered and filtered-transposed matrix for each input file  

---

### Baseline Command (Example Run)

#### (1) Load Java on HPC cluster
```bash
module load ccs/java/java-11.0.2
```

#### (2) Compile the Java program
```bash
javac ExpressionMatrixProcessor.java
```

#### (3) Run on gene-level data
```bash
nohup java ExpressionMatrixProcessor --input-dir RawData > ExpressionMatrixProcessor_gene.log 2>&1 &
```

#### (4) (Optional) Monitor running process
```bash
ps aux | grep ExpressionMatrixProcessor
tail -f ExpressionMatrixProcessor_gene.log
```

#### (5) Run on isoform-level data
```bash
nohup java ExpressionMatrixProcessor --two-column-header --input-dir RawData_iso > ExpressionMatrixProcessor_iso.log 2>&1 &
```

#### (6) (Optional) Monitor running process
```bash
ps aux | grep ExpressionMatrixProcessor
tail -f ExpressionMatrixProcessor_iso.log
```

---

### Input Requirements

- Input directory must contain one or more `.txt` raw count matrices.  
- Only files ending with `.txt` are processed; others are ignored.  
- Include only `.txt` files you wish to process.  
- Isoform and gene files must be placed in **separate directories** (e.g., `RawData_gene/` and `RawData_iso/`).  
- Matrices are typically large and should be processed on an **HPC cluster**.

#### Header Formats

**(1) One-column header (gene-level)**  
```
GeneID    Cell1    Cell2 ...
ENSGXX    0        1.333333 ...
```

**(2) Two-column header (isoform-level)**  
```
TranscriptID  GeneID  Cell1  Cell2 ...
ENSTXX        ENSGXX  0      2.33333 ...
```

Use `--two-column-header` for isoform-level input files.

---

### Command-Line Options

| Option | Description | Default |
|--------|--------------|----------|
| `--two-column-header` | Indicates input has `TranscriptID` and `GeneID` columns. Used for isoform-level data. | Off |
| `--input-dir <DIR>` | Directory containing `.txt` matrices. | `RawData` |
| `--cell-threshold <DOUBLE>` | Minimum total reads per cell to retain. | 500 |
| `--gene-min-cells <INT>` | Minimum number of cells where a gene/isoform must be expressed (>0) to retain. | 10 |

**Examples**
```bash
java ExpressionMatrixProcessor --input-dir RawData_gene --cell-threshold 300 --gene-min-cells 10 > gene.log 2>&1
java ExpressionMatrixProcessor --two-column-header --input-dir RawData_iso --cell-threshold 300 --gene-min-cells 10 > iso.log 2>&1
```

---

### 📤 Output Files (per input `Sample1.txt`)

| Output File | Description |
|--------------|-------------|
| `Sample1.filtered_expression_matrix.txt` | Filtered file after removing low-depth cells. |
| `Sample1.filtered_transposed_expression_matrix.txt` | Transposed version (rows = CellID, columns = genes/isoforms). |
| `removed_cells_Sample1.txt` | List of cells removed due to low total read counts. |
| `filtered_genes.txt` | List of genes/isoforms removed due to low expression (once per run). |
| `Sample1.null_values.txt` | Debug output file (normally empty). |

**Temporary Files**  
- Created under `temp_transpose_<sample>/` during transposition.

---

### Acknowledgment

Developed by Mark Ebbert, University of Kentucky.  
Intended for use in long-read single-cell and isoform-level RNA expression processing pipelines.  
If used in a publication, please cite the associated paper, methods, and repository.
