package main

// ... existing imports ...
import (
	// Add this import
	"fmt"
	"log"
	"os"
)

// ... existing functions ...

func main() {
	if len(os.Args) != 3 {
		fmt.Println("Usage: ./preprocess <dataset_type> <file_path>")
		fmt.Println("dataset_type: 'rat' or 'golub'")
		os.Exit(1)
	}

	datasetType := os.Args[1]
	filePath := os.Args[2]

	switch datasetType {
	case "rat":
		// Process rat data (GDS2901.soft)
		dataWithGenes, err := ReadData(filePath)
		if err != nil {
			log.Fatalf("Error reading rat data: %v", err)
		}

		// Remove last row and probeset 2475
		dataWithGenes.Data = removeRow(dataWithGenes.Data, len(dataWithGenes.GeneIDs)-1)
		dataWithGenes.GeneIDs = append(dataWithGenes.GeneIDs[:2474], dataWithGenes.GeneIDs[2475:]...)
		dataWithGenes.Data = removeRow(dataWithGenes.Data, 2474)

		// Process data
		logData := applyLog2(dataWithGenes.Data)
		normData := NormalizeQuantiles(logData)

		// Extract conditions
		ekerMutants := ExtractEkerMutants(normData)
		wildTypes := extractWildTypes(normData)

		// Save with gene IDs
		if err := saveToCSV(ekerMutants, dataWithGenes.GeneIDs, "eker_mutants.csv"); err != nil {
			log.Fatalf("Error saving Eker mutants data: %v", err)
		}
		if err := saveToCSV(wildTypes, dataWithGenes.GeneIDs, "wild_types.csv"); err != nil {
			log.Fatalf("Error saving wild types data: %v", err)
		}
		fmt.Println("Rat data processing complete! Files saved: eker_mutants.csv, wild_types.csv")

	case "golub":
		// Process Golub data
		dataWithGenes, err := ReadGolubData(filePath)
		if err != nil {
			log.Fatalf("Error reading Golub data: %v", err)
		}

		// Process data
		logData := applyLog2(dataWithGenes.Data)
		normData := NormalizeQuantiles(logData)

		// Extract ALL and AML samples (using same functions as Eker/Wild type)
		allSamples := ExtractALLSamples(normData) // First condition (ALL)
		amlSamples := ExtractAMLSamples(normData) // Second condition (AML)

		// Save with gene IDs
		if err := saveToCSV(allSamples, dataWithGenes.GeneIDs, "all_samples.csv"); err != nil {
			log.Fatalf("Error saving ALL samples data: %v", err)
		}
		if err := saveToCSV(amlSamples, dataWithGenes.GeneIDs, "aml_samples.csv"); err != nil {
			log.Fatalf("Error saving AML samples data: %v", err)
		}
		fmt.Println("Golub data processing complete! Files saved: all_samples.csv, aml_samples.csv")

	default:
		log.Fatalf("Unknown dataset type: %s. Use 'rat' or 'golub'", datasetType)
	}
}
