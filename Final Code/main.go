package main

import (
	"fmt"
	"log"
	"os"

	"gonum.org/v1/gonum/mat"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Println("Usage: ./preprocess <dataset_type> <file_path>")
		fmt.Println("dataset_type: 'rat' or 'golub'")
		os.Exit(1)
	}

	datasetType := os.Args[1]
	filePath := os.Args[2]

	// Create output directories if they don't exist
	for _, dir := range []string{"output/diffcoex", "output/coxpress"} {
		if err := os.MkdirAll(dir, 0755); err != nil {
			log.Fatalf("Error creating output directory %s: %v", dir, err)
		}
	}

	// Process data using both methods
	switch datasetType {
	case "rat":
		if err := processRatData(filePath); err != nil {
			log.Fatalf("Error processing rat data: %v", err)
		}
		fmt.Println("Rat data processing complete! Files saved:")
		fmt.Println("- output/diffcoex/rat_eker_mutants.csv")
		fmt.Println("- output/diffcoex/rat_wild_types.csv")
		fmt.Println("- output/coxpress/rat_eker_mutants.csv")
		fmt.Println("- output/coxpress/rat_wild_types.csv")

	case "golub":
		if err := processGolubData(filePath); err != nil {
			log.Fatalf("Error processing Golub data: %v", err)
		}
		fmt.Println("Golub data processing complete! Files saved:")
		fmt.Println("- output/diffcoex/golub_ALL_samples.csv")
		fmt.Println("- output/diffcoex/golub_AML_samples.csv")
		fmt.Println("- output/coxpress/golub_ALL_samples.csv")
		fmt.Println("- output/coxpress/golub_AML_samples.csv")

	default:
		log.Fatalf("Unknown dataset type: %s. Use 'rat' or 'golub'", datasetType)
	}
}

func processRatData(filePath string) error {
	fmt.Printf("Reading rat data from: %s\n", filePath)

	// Read data
	dataWithGenes, err := ReadData(filePath)
	if err != nil {
		return fmt.Errorf("error reading rat data: %v", err)
	}

	// DiffCoEx preprocessing
	{
		// Create a deep copy for DiffCoEx processing
		diffCoExData := copyDataWithGenes(dataWithGenes)

		// Remove last row and probeset 2475
		diffCoExData.Data = removeRow(diffCoExData.Data, len(diffCoExData.GeneIDs)-1)
		diffCoExData.GeneIDs = append(diffCoExData.GeneIDs[:2474], diffCoExData.GeneIDs[2475:]...)
		diffCoExData.Data = removeRow(diffCoExData.Data, 2474)

		// Process data
		logData := applyLog2(diffCoExData.Data)
		normData := NormalizeQuantiles(logData)

		// Extract conditions
		ekerMutants := ExtractEkerSamples(normData)
		wildTypes := ExtractWildSamples(normData)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(ekerMutants, diffCoExData.GeneIDs, "output/diffcoex/rat_eker_mutants.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx Eker mutants: %v", err)
		}
		if err := saveToCSV(wildTypes, diffCoExData.GeneIDs, "output/diffcoex/rat_wild_types.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx wild types: %v", err)
		}
	}

	// coXpress preprocessing
	{
		// Extract conditions directly from raw data
		ekerMutants := ExtractEkerSamples(dataWithGenes.Data)
		wildTypes := ExtractWildSamples(dataWithGenes.Data)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(ekerMutants, dataWithGenes.GeneIDs, "output/coxpress/rat_eker_mutants.csv"); err != nil {
			return fmt.Errorf("error saving coXpress Eker mutants: %v", err)
		}
		if err := saveToCSV(wildTypes, dataWithGenes.GeneIDs, "output/coxpress/rat_wild_types.csv"); err != nil {
			return fmt.Errorf("error saving coXpress wild types: %v", err)
		}
	}

	return nil
}

func processGolubData(filePath string) error {
	fmt.Printf("Reading Golub data from: %s\n", filePath)

	// Read data
	dataWithGenes, err := ReadGolubData(filePath)
	if err != nil {
		return fmt.Errorf("error reading Golub data: %v", err)
	}

	// DiffCoEx preprocessing
	{
		// Create a deep copy for DiffCoEx processing
		diffCoExData := copyDataWithGenes(dataWithGenes)

		// Process data
		logData := applyLog2(diffCoExData.Data)
		normData := NormalizeQuantiles(logData)

		// Extract ALL and AML samples
		allSamples := ExtractALLSamples(normData)
		amlSamples := ExtractAMLSamples(normData)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(allSamples, diffCoExData.GeneIDs, "output/diffcoex/golub_ALL_samples.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx ALL samples: %v", err)
		}
		if err := saveToCSV(amlSamples, diffCoExData.GeneIDs, "output/diffcoex/golub_AML_samples.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx AML samples: %v", err)
		}
	}

	// coXpress preprocessing
	{
		// Extract ALL and AML samples directly from raw data
		allSamples := ExtractALLSamples(dataWithGenes.Data)
		amlSamples := ExtractAMLSamples(dataWithGenes.Data)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(allSamples, dataWithGenes.GeneIDs, "output/coxpress/golub_ALL_samples.csv"); err != nil {
			return fmt.Errorf("error saving coXpress ALL samples: %v", err)
		}
		if err := saveToCSV(amlSamples, dataWithGenes.GeneIDs, "output/coxpress/golub_AML_samples.csv"); err != nil {
			return fmt.Errorf("error saving coXpress AML samples: %v", err)
		}
	}

	return nil
}

func copyDataWithGenes(d *DataWithGenes) *DataWithGenes {
	newData := mat.NewDense(d.Data.RawMatrix().Rows, d.Data.RawMatrix().Cols, nil)
	newData.Copy(d.Data)

	newGeneIDs := make([]string, len(d.GeneIDs))
	copy(newGeneIDs, d.GeneIDs)

	return &DataWithGenes{
		Data:    newData,
		GeneIDs: newGeneIDs,
	}
}
