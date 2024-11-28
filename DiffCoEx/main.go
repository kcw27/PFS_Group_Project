package main

import (
	"fmt"
	"log"
)

func main() {
	dataFile := "GDS2901.soft"

	rawData, err := ReadData(dataFile)
	if err != nil {
		log.Fatalf("Error reading data: %v", err)
	}

	// Remove probeset 2475 (as done in R code)
	rawDataWithoutProbe := removeRow(rawData, 2474)

	// Log2 transformation
	logData := applyLog2(rawDataWithoutProbe)

	// Quantile normalization
	normData := NormalizeQuantiles(logData)

	// Extract Eker mutants and wild-types
	ekerMutants := ExtractEkerMutants(normData)
	wildTypes := extractWildTypes(normData)

	// Save Eker mutants to CSV
	if err := saveToCSV(ekerMutants, "eker_mutants.csv"); err != nil {
		log.Fatalf("Error saving Eker mutants data: %v", err)
	}

	// Save wild types to CSV
	if err := saveToCSV(wildTypes, "wild_types.csv"); err != nil {
		log.Fatalf("Error saving wild types data: %v", err)
	}

	fmt.Println("Data processing complete. Files saved: eker_mutants.csv, wild_types.csv")

	// NEED actual colorhC1C2, implementation below is example data
	colorh1C1C2 := map[string]string{"gene1": "red", "gene2": "blue"}

	// Generate permutations
	numPermutations := 1000
	permutations := generatePermutations(ekerMutants, wildTypes, numPermutations)

	// Scale and combine the data
	d := combineAndScaleData(ekerMutants, wildTypes)

	// Compute dispersion matrix and null distribution
	uniqueColors := []string{"red", "blue"} // Example unique colors
	dispersionMatrix := make([][]float64, len(uniqueColors))
	nullDistrib := make(map[string]map[string][]float64)

	for i, c1 := range uniqueColors {
		dispersionMatrix[i] = make([]float64, len(uniqueColors))
		nullDistrib[c1] = make(map[string][]float64)
		for j, c2 := range uniqueColors {
			// Calculate dispersion
			dispersionMatrix[i][j] = dispersionModule2Module(c1, c2, ekerMutants, wildTypes, colorh1C1C2)

			// Generate null distribution
			nullDistrib[c1][c2] = make([]float64, numPermutations)
			for k, perm := range permutations {
				nullDistrib[c1][c2][k] = permutationProcedureModule2Module(perm, d, c1, c2, colorh1C1C2)
			}
		}
	}

	// Create permutation summary
	permutationSummary := make([][]int, len(uniqueColors))
	for i := range permutationSummary {
		permutationSummary[i] = make([]int, len(uniqueColors))
		for j := range permutationSummary[i] {
			count := 0
			for _, val := range nullDistrib[uniqueColors[i]][uniqueColors[j]] {
				if val >= dispersionMatrix[i][j] {
					count++
				}
			}
			permutationSummary[i][j] = count
		}
	}

	// Print results
	fmt.Println("Dispersion Matrix:", dispersionMatrix)
	fmt.Println("Permutation Summary:", permutationSummary)

}
