package main

import (
	"fmt"
	"log"
	"os"
)

func main() {
	// These would come from R Shiny as command-line arguments
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
		ekerMutants, wildTypes, err := preprocessForCoXpress(filePath)
		if err != nil {
			log.Fatalf("Error preprocessing rat data: %v", err)
		}

		// Save the processed data to files
		if err := saveToCSV(ekerMutants, "eker_mutants.csv"); err != nil {
			log.Fatalf("Error saving Eker mutants data: %v", err)
		}
		if err := saveToCSV(wildTypes, "wild_types.csv"); err != nil {
			log.Fatalf("Error saving wild types data: %v", err)
		}
		fmt.Println("Rat data processing complete! Files saved: eker_mutants.csv, wild_types.csv")

	case "golub":
		// Process Golub data
		rawData, err := ReadGolubData(filePath)
		if err != nil {
			log.Fatalf("Error reading Golub data: %v", err)
		}

		// Extract ALL and AML samples
		allSamples := ExtractALLSamples(rawData)
		amlSamples := ExtractAMLSamples(rawData)

		// Save ALL samples to CSV
		if err := saveToCSV(allSamples, "all_samples.csv"); err != nil {
			log.Fatalf("Error saving ALL samples data: %v", err)
		}
		// Save AML samples to CSV
		if err := saveToCSV(amlSamples, "aml_samples.csv"); err != nil {
			log.Fatalf("Error saving AML samples data: %v", err)
		}
		fmt.Println("Golub data processing complete! Files saved: all_samples.csv, aml_samples.csv")

	default:
		log.Fatalf("Unknown dataset type: %s. Use 'rat' or 'golub'", datasetType)
	}
}
