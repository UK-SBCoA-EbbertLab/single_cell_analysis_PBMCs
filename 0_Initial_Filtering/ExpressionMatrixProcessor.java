// -*- coding: utf-8 -*-
// package ebbertLab.expressionMatrixProcessor;

import java.io.*;
import java.util.*;

public class ExpressionMatrixProcessor {

    // Configuration constants for controlling filtering thresholds and input directory
//    static boolean twoColumnHeader = false;       // Declare whether input files have one or two header columns (i.e., gene ID, or gene and isoform ID)
//    static final double CELL_THRESHOLD = 500.0;   // Minimum sum of reads for retaining a cell (column filtering threshold)
//    static final int GENE_MIN_CELLS = 10;         // Minimum number of cells in which a gene must be expressed (>0) to be retained
//    static final String INPUT_DIR = "RawData";    // Directory path where input expression matrix files (.txt format) are stored

    public static void main(String[] args) throws IOException {

    	// === Command-line Argument Parsing ===
    	boolean twoColumnHeader = false;    // Default
    	String inputDir = "RawData";        // Default
    	double cellThreshold = 500.0;       // Default
    	int geneMinCells = 10;              // Default

    	for (int i = 0; i < args.length; i++) {
    	    switch (args[i]) {
    	        case "--two-column-header":
    	            twoColumnHeader = true;
    	            break;

    	        case "--input-dir":
    	            if (i + 1 < args.length) {
    	                inputDir = args[++i];
    	            } else {
    	                System.out.println("❌ Missing value for --input-dir");
    	                System.exit(1);
    	            }
    	            break;

    	        case "--cell-threshold":
    	            if (i + 1 < args.length) {
    	                cellThreshold = Double.parseDouble(args[++i]);
    	            } else {
    	                System.out.println("❌ Missing value for --cell-threshold");
    	                System.exit(1);
    	            }
    	            break;

    	        case "--gene-min-cells":
    	            if (i + 1 < args.length) {
    	                geneMinCells = Integer.parseInt(args[++i]);
    	            } else {
    	                System.out.println("❌ Missing value for --gene-min-cells");
    	                System.exit(1);
    	            }
    	            break;

    	        default:
    	            System.out.println("⚠️ Unknown argument: " + args[i]);
    	            break;
    	    }
    	}

        System.out.println("🔍 Header mode: " + (twoColumnHeader ? "TWO columns (gene name + ID)" : "ONE column (gene name only)"));
            
		/*
		 * If two-column header file, the expression column starts at index 2. Otherwise, 1.
		 */
		int expressionColumnIndex;
		if (twoColumnHeader) {
			expressionColumnIndex = 2;
		} else {
			expressionColumnIndex = 1;
		}

        System.out.println("🔍 Searching for .txt expression matrices in directory: " + inputDir);

        // Step 1: Locate input files
        File dir = new File(inputDir);
        File[] inputFiles = dir.listFiles((d, name) -> name.endsWith(".txt"));

        if (inputFiles == null || inputFiles.length == 0) {
            System.out.println("❌ No .txt files found in directory: " + inputDir);
            return; // Terminate if no input files are found
        }

        System.out.println("📂 Found " + inputFiles.length + " input files to process.");

        // === Data Structures ===
        Map<String, Integer> globalGeneCounts = new LinkedHashMap<>();  // Track total cells expressing each gene
        Map<File, String[]> geneOrderPerFile = new HashMap<>();         // Remember gene order for each file
        Map<File, File> filteredFilePerInput = new HashMap<>();         // Map each input file to its filtered output file

        long pipelineStart = System.nanoTime();  // Timer start

        // ============================ PHASE 1 ============================
        // Cell filtering and global gene expression counting phase
        System.out.println("\n=== PHASE 1: FILTER CELLS + TRACK GENE COUNTS ===");

        // Process each input file sequentially
        for (File inputFile : inputFiles) {
            String sampleName = inputFile.getName().replace(".txt", "");
            System.out.println("\n🔄 Starting file: " + sampleName);
            long fileStart = System.nanoTime();

            // Track filtered-out cells for this file
            BufferedWriter filteredCellsWriter = new BufferedWriter(new FileWriter("removed_cells_" + sampleName + ".txt"));
            filteredCellsWriter.write("CellID\tTotalReads\n");

            // First Pass: Compute total counts for each cell (sum columns)
            BufferedReader reader = new BufferedReader(new FileReader(inputFile));
            String headerLine = reader.readLine(); // Read header line to extract column identifiers
            String[] header = headerLine.split("\t");
            
            int headerColumns = twoColumnHeader ? 2 : 1;     // <-- key branching logic
            int numColumns = header.length - headerColumns;  // <-- adjusts based on header type
//            int numColumns = header.length - 1;    // Deduct gene 


			
            /*
             * Get cellIDs. If two-column header file, start at the third column (index 2)
             */
			String[] cellIDs = Arrays.copyOfRange(header, expressionColumnIndex, header.length); // Cell IDs
			
            double[] columnSums = new double[numColumns]; // Initialize sums array
            Arrays.fill(columnSums, 0.0);

            System.out.println("  🔬 First pass: Computing column sums for " + numColumns + " cells...");

            // Stream through the file to compute column sums incrementally
            String line;
            int geneRows = 0;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                for (int i = 0; i < numColumns; i++) {
                    columnSums[i] += Double.parseDouble(fields[i + expressionColumnIndex]);
                }
                geneRows++;
                if (geneRows % 5000 == 0) {
                    System.out.println("    Processed " + geneRows + " gene rows...");
                }
            }
            reader.close();
            System.out.println("  ✅ Finished computing column sums.");

            // Determine which cells to retain based on cellThreshold
            List<Integer> retainedIndices = new ArrayList<>();
			System.out.println("  ❌ Filtering low-count cells (< " + cellThreshold + " reads):");
            int cellsFiltered = 0;
            for (int i = 0; i < numColumns; i++) {
                if (columnSums[i] >= cellThreshold) {
                    retainedIndices.add(i);
                } else {
                    filteredCellsWriter.write(cellIDs[i] + "\t" + columnSums[i] + "\n");
                    cellsFiltered++;
                }
            }
            filteredCellsWriter.close();
            System.out.println("  📊 Cells retained: " + retainedIndices.size() + " | Filtered out: " + cellsFiltered);

            
            
            System.out.println("  🛠️ Second pass: Writing filtered matrix + counting gene expression...");
            
            // Second Pass: Write filtered matrix + update global gene counts
            File filteredOutput = new File(sampleName + ".filtered_expression_matrix.txt");
            BufferedWriter writer = new BufferedWriter(new FileWriter(filteredOutput));

//            writer.write(header[0]);  // Write gene name 
            if (twoColumnHeader) {
                writer.write(header[0] + "\t" + header[1]);
            } else {
                writer.write(header[0]);
            }
            for (int idx : retainedIndices) {
                writer.write("\t" + header[idx + expressionColumnIndex]);
            }
            writer.newLine();

            List<String> geneOrder = new ArrayList<>();
            int genesProcessed = 0;

            // Begin file again and advance past header row
            reader = new BufferedReader(new FileReader(inputFile));
            headerLine = reader.readLine();

            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
//                String geneName = fields[0];
                
                /*
                 * If a two-column header file, use both gene name and tx ID as geneName
                 */
                String geneName = twoColumnHeader ? fields[0] + "|" + fields[1] : fields[0];
                geneOrder.add(geneName);

                int expressedInCells = 0;

//                writer.write(fields[0]);
                if (twoColumnHeader) {
                    writer.write(fields[0] + "\t" + fields[1]);
                } else {
                    writer.write(fields[0]);
                }

                /*
                 * Print all of the expression values
                 */
                for (int idx : retainedIndices) {
                	
                	/*
                	 * Round to two decimals before assessing if it's a whole
                	 * number and printing. A lot of numbers are N.000001. I want
                	 * those to round to the whole number and then be printed as
                	 * a whole number in the next line of code.
                	 */

					double value = Math.round(Double.parseDouble(fields[idx + expressionColumnIndex]) * 100.0) / 100.0;
                    // writer.write("\t" + value);
					
					// If it's a whole number, don't print the .0 to save space
					writer.write("\t" + (value == (long) value ? String.valueOf((long) value) : String.format("%.2f", value)));
					
                    if (value > 0) {
                        expressedInCells++;
                    }
                }
                writer.newLine();

                globalGeneCounts.put(geneName, globalGeneCounts.getOrDefault(geneName, 0) + expressedInCells);

                genesProcessed++;
                if (genesProcessed % 5000 == 0) {
                    System.out.println("    Processed " + genesProcessed + " genes...");
                }
            }

            reader.close();
            writer.close();

            System.out.println("  ✅ Filtered matrix written: " + filteredOutput.getName());

            // Record metadata for downstream use
            geneOrderPerFile.put(inputFile, geneOrder.toArray(new String[0]));
            filteredFilePerInput.put(inputFile, filteredOutput);

            long fileEnd = System.nanoTime();
            System.out.printf("⏱️ Completed file: %s in %.2f seconds\n", sampleName, (fileEnd - fileStart) / 1_000_000_000.0);
        }
        
        
        /*
         * Verify gene order is identical across all input files.
         */
        if(areAllArraysIdentical(geneOrderPerFile)){
        	System.out.println("\n✅ All files have identical genes and order.");
        }
        else {
        	System.out.println("❌ ERROR: Some files had different genes / order.");
        }
        

        // ======================== PHASE 2 ========================
        // Global gene filtering based on expression frequency across all input files
        System.out.println("\n=== PHASE 2: Global gene filtering based on expression counts ===");
		long geneFilterStart = System.nanoTime();

        BufferedWriter filteredGenesWriter = new BufferedWriter(new FileWriter("filtered_genes.txt"));
        filteredGenesWriter.write("GeneID\tExpressingCells\n");
        Set<String> globallyFilteredGenes = new HashSet<>();

        // Determine genes to be globally filtered and write them to output
        int genesFiltered = 0;
        for (Map.Entry<String, Integer> entry : globalGeneCounts.entrySet()) {
            String gene = entry.getKey();
            int count = entry.getValue();
            if (count < geneMinCells) {
                globallyFilteredGenes.add(gene);
                filteredGenesWriter.write(gene + "\t" + count + "\n");
                // System.out.println("  🔥 Filtering out gene: " + gene + " (Expressed in " + count + " cells)");
                genesFiltered++;
            }
        }
        filteredGenesWriter.close();
        System.out.println("  🚫 Total globally filtered genes: " + genesFiltered);
		long geneFilterEnd = System.nanoTime();
		System.out.printf("⏱️ Completed global gene filtering: in %.2f seconds\n", (geneFilterEnd - geneFilterStart) / 1_000_000_000.0);

		
		
		
        // ======================== PHASE 3 ========================
        // Transpose filtered matrices while omitting globally filtered genes
        System.out.println("\n=== PHASE 3: Transposing filtered matrices ===");

        /*
         * For each original input file:
         * 	1. Get the filtered file and gene order
         * 	2. Print each remaining cell to separate files while omitting genes that
         *     where globally omitted.
         */
        for (File inputFile : inputFiles) {
            String sampleName = inputFile.getName().replace(".txt", "");
            System.out.println("\n🔄 Transposing filtered matrix for file: " + sampleName);
            long transposeStart = System.nanoTime();

            // Get file names for filtered file and gene order for each original input file
            File filteredFile = filteredFilePerInput.get(inputFile);
            String[] geneOrder = geneOrderPerFile.get(inputFile);
            
            // There should be geneOrder.length - globallyFilteredGenes.size() genes printed,
            // in the end. ChatGPT got this wrong originally. Was printing only oeneOrder.length
            System.out.println("\n  There are " + (geneOrder.length - globallyFilteredGenes.size()) + " genes remaining.");

            BufferedReader reader = new BufferedReader(new FileReader(filteredFile));
            String headerLine = reader.readLine();
            String[] header = headerLine.split("\t");
            String[] cellIDs = Arrays.copyOfRange(header, expressionColumnIndex, header.length);

            
            System.out.println("\n  Writing each remaining cell to unique tmp file, while ignoring globally eliminated genes.");
            /* 
             * Create writers for temporary files to build transposed matrix. Write each cell to 
             * a unique file that can be accessed when transposing.
             */
            String tempDir = "temp_transpose_" + sampleName;
            new File(tempDir).mkdir();
            BufferedWriter[] tempWriters = new BufferedWriter[cellIDs.length];
            for (int i = 0; i < cellIDs.length; i++) {
                tempWriters[i] = new BufferedWriter(new FileWriter(tempDir + "/" + cellIDs[i] + ".tmp"));
            }


            /* 
             * For current original input file, loop over its filtered file
             * (the file excluding cells with <500 reads)
             * and write each cell to its own file while skipping filtered genes
             * (genes that were not expressed in at least 10 cells across all
             * files).
             */
            int genesWritten = 0;
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
//                String gene = fields[0]; // Get gene name
                String gene = twoColumnHeader ? fields[0] + "|" + fields[1] : fields[0]; // Get gene name
                if (globallyFilteredGenes.contains(gene)) {
                    continue;  // Skip globally filtered gene
                }
                for (int i = 0; i < cellIDs.length; i++) {
                    tempWriters[i].write(fields[i + expressionColumnIndex]);
                    tempWriters[i].newLine();
                }
                genesWritten++;
                if (genesWritten % 5000 == 0) {
                    System.out.println("  Writing gene row: " + genesWritten);
                }
            }
            reader.close();
            for (BufferedWriter tempWriter : tempWriters) {
                tempWriter.close();
            }

            /*
             * Compose final transposed matrix from temporary files
             */
            System.out.println("  📊 Writing final transposed matrix...");

            BufferedWriter writer = new BufferedWriter(new FileWriter(sampleName + ".filtered_transposed_expression_matrix.txt"));
            BufferedWriter debugWriter = new BufferedWriter(new FileWriter(sampleName + ".null_values.txt"));
            
            /*
             * Write header row, which is now the gene names
             */
            writer.write("CellID");
            for (String gene : geneOrder) {
            	/*
            	 * Ignore globallyFilteredGenes
            	 */
                if (!globallyFilteredGenes.contains(gene)) {
                    writer.write("\t" + gene);
                }
            }
            writer.newLine();

            /*
             * For each cell, write the values for all genes, ignoring those that
             * were globally filtered.
             */
            for (int i = 0; i < cellIDs.length; i++) {
                writer.write(cellIDs[i]); // Write the cell ID (first item in row, which will form the header column)
                
                /*
                 * Open the temp file containing gene expression values for this cell
                 */
                BufferedReader tempReader = new BufferedReader(new FileReader(tempDir + "/" + cellIDs[i] + ".tmp"));
                int geneCount = 0;
                String value;
                
                /*
                 * For each gene found in the original file, print gene expression
                 * values for genes that were not globally excluded.
                 */
                for (String gene : geneOrder) {
                	geneCount++;
                    if (globallyFilteredGenes.contains(gene)) {
                    	
                    	/*
                    	 * ChatGPT originally advanced the reader to the next line
                    	 * if the gene was globally filtered, but the files we're 
                    	 * reading from already do not have these values in them
                    	 * because the globally filtered. This should be why we
                    	 * were getting null values. i.e., the reader gets out
                    	 * of sync by advancing when it does not. readLine() returns
                    	 * null when it runs out of text.
                    	 */
//                        tempReader.readLine(); 
                        continue;
                    }
                    value = tempReader.readLine();
                    
                    if (value == null) {
                    	System.out.println("ERROR: tmpeReader.readLine() returned null."); 
                    	System.exit(-1);
                    }
                    writer.write("\t" + value);
                    
                    if(null == value) {
                    	debugWriter.write("Gene count: " + geneCount + "; Gene: " + gene + "; Value: " + value + "\n");
                    }
                }
                tempReader.close();
                writer.newLine();
                if ((i + 1) % 5000 == 0) {
                    System.out.println("  Completed " + (i + 1) + " of " + cellIDs.length + " transposed rows.");
                }
            }
            writer.close();
            debugWriter.close();

            // Clean up temporary files
            for (File file : new File(tempDir).listFiles()) {
//                file.delete();
            }
//            new File(tempDir).delete();

            long transposeEnd = System.nanoTime();
            System.out.printf("  ✅ Transposition complete for %s in %.2f seconds\n", sampleName, (transposeEnd - transposeStart) / 1_000_000_000.0);
        }

        // Pipeline complete
        long pipelineEnd = System.nanoTime();
        System.out.printf("\n🎉 Pipeline complete for all files in %.2f seconds\n", (pipelineEnd - pipelineStart) / 1_000_000_000.0);
    }
    

    
    /*
     * Verify all arrays are identical
     */
    public static boolean areAllArraysIdentical(Map<File, String[]> map) {

    	/*
    	 * If no arrays provided, throw error
    	 */
    	if (map == null || map.isEmpty()) {
    		System.out.println("No arrays provided!");
    		return false;
    	}

    	File referenceFile = map.keySet().iterator().next();
    	String[] referenceArray = map.get(referenceFile);

    	for (Map.Entry<File, String[]> entry : map.entrySet()) {
    		File currentFile = entry.getKey();
    		String[] currentArray = entry.getValue();

    		if (referenceArray.length != currentArray.length) {
    			System.out.println("Mismatch in length. \n" +
    					referenceFile.getName() + " length: " + referenceArray.length + "\n" +
    					currentFile.getName() + " length: " + currentArray.length);
    			return false;
    		}

    		for (int i = 0; i < referenceArray.length; i++) {
    			if (!referenceArray[i].equals(currentArray[i])) {
    				System.out.println("Mismatch in file: " + currentFile.getName() + " at index " + i);

    				System.out.println("Mismatch in order or gene name at index " + i + ". \n" +
    						referenceFile.getName() + " gene: " + referenceArray[i] + "\n" +
    						currentFile.getName() + " gene: " + currentArray[i]);
    				return false;
    			}
    		}
    	}
    	return true;
    }
}
