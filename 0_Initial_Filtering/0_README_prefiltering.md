## ExpressionMatrixProcessor
This Java processing script is designed to filter and transpose tab-delimited expression matrices (genes/isoforms x cells) to reduce file size for downstream supervised filtering. It performs the following functions: 
(1) Removes low-depth cells
(2) Removes rarely expressed genes and isoforms
(3) Rounds numeric values to two decimal places. If the result after rounding is an integer, decimals are removed.
(4) Write a filtered and filtered-transposed matrix for each input file. 

### BASELINE COMMAND (example run)
# (1) load java on HPC cluster
module load ccs/java/java-11.0.2

# (2) Compile and run on HPC node
javac ExpressionMatrixProcessor.java

# (3) Run on gene-level data
nohup java ExpressionMatrixProcessor --input-dir RawData > ExpressionMatrixProcessor_gene.log 2>&1 &

# (4) (optional) Monitor running process
ps aux | grep ExpressionMatrixProcessor
tail -f ExpressionMatrixProcessor_gene.log

# (5) Run on isoform-level data
nohup java ExpressionMatrixProcessor --two-column-header --input-dir RawData_iso > ExpressionMatrixProcessor_iso.log 2>&1 &

# (6) (optional) Monitor running process
ps aux | grep ExpressionMatrixProcessor
tail -f ExpressionMatrixProcessor_iso.log


### INPUT REQUIREMENTS
- Input directory must contain at least one .txt raw count matrices.
- The program only processes files ending with .txt, all others are ignored
- Only include .txt files in the input directory that you want processed. 
- Isoform and gene files must be separated in distinct directories (e.g., RawData_gene and RawData_iso) within the home one, since each will require a separate iteration. 
- These input matrices can be very large and will typically require execution on a high-performance computing (HPC) cluster.

Header formats:
(1) One-column header (gene-level)
    GeneID    Cell1    Cell2 ...
    ENSGXX    0        1.333333  ...

(2) Two-column header (isoform-level)
    TranscriptID  GeneID  Cell1  Cell2 ...
    ENSTXX        ENSGXX  0      2.33333  ...

Use --two-column-header for isoform-level input files.


### COMMAND-LINE OPTIONS
--two-column-header
    Indicates input has TranscriptID and GeneID columns (no argument). This is used for isoform-level data and is not needed to process gene-level data. (Default = Off)

--input-dir <DIR>
    Directory containing .txt matrices (default = RawData)

--cell-threshold <DOUBLE>
    Minimum total reads per cell to retain (default = 500)

--gene-min-cells <INT>
    Minimum number of cells where gene/isoform must be expressed (expression count >0) to be retained (default = 10)

Examples:
java ExpressionMatrixProcessor --input-dir RawData_gene --cell-threshold 300 --gene-min-cells 10 > gene.log 2>&1
java ExpressionMatrixProcessor --two-column-header --input-dir RawData_iso --cell-threshold 300 --gene-min-cells 10 > iso.log 2>&1


### OUTPUT FILES (per input Sample1.txt)
Sample1.filtered_expression_matrix.txt
    Filtered file after removing low-depth cells.

Sample1.filtered_transposed_expression_matrix.txt
    Transposed version (rows = CellID, columns = genes/isoforms).

removed_cells_Sample1.txt
    List of cells removed due to low total read counts.

filtered_genes.txt
    List of genes/isoforms removed due to low expression (once per run).

Sample1.null_values.txt
    Debug output file (should normally be empty).

Temporary files:
    Created under temp_transpose_<sample>/ during transposition.


#### ACKNOWLEDGMENT
Developed by Mark Ebbert, University of Kentucky. 
Intended for use in long-read single-cell and isoform-level RNA expression processing pipelines. If used in a publication, please cite the associated paper, methods, and repository.
