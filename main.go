package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strconv"
	"math"
	"gonum.org/v1/gonum/stat"
)

type ModuleStats struct {
	Name string
	TStatistic float64
	PValue float64
	Size int
}

func main() {
	// Load module assignments
	moduleMap, err := loadModules("data/golub/golub_diffcoex.csv")
	if err != nil {
		log.Fatal("Error loading modules:", err)
	}

	// Load expression data
	amlData, err := loadExpressionData("data/golub/aml_samples.csv")
	if err != nil {
		log.Fatal("Error loading AML data:", err)
	}

	allData, err := loadExpressionData("data/golub/all_samples.csv")
	if err != nil {
		log.Fatal("Error loading ALL data:", err)
	}

	// Analyze each module
	fmt.Printf("Module\tSize\tT-Statistic\tP-Value\n")
	for module := range getUniqueModules(moduleMap) {
		stats := analyzeModule(module, moduleMap, amlData, allData)
		fmt.Printf("%s\t%d\t%f\t%f\n", stats.Name, stats.Size, stats.TStatistic, stats.PValue)
	}
}

func loadModules(filename string) (map[string]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	// Skip header
	_, err = reader.Read()
	if err != nil {
		return nil, err
	}

	moduleMap := make(map[string]string) // gene -> module
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		moduleMap[record[0]] = record[1]
	}
	return moduleMap, nil
}

func loadExpressionData(filename string) (map[string][]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	data := make(map[string][]float64)
	
	// Skip header if it exists
	_, err = reader.Read()
	if err != nil {
		return nil, err
	}
	
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		
		geneName := record[0]
		values := make([]float64, 0)
		
		// Convert string values to float64, starting from column 1
		for _, val := range record[1:] {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				continue
			}
			values = append(values, f)
		}
		
		if len(values) > 0 {
			data[geneName] = values
		}
	}
	return data, nil
}

func analyzeModule(moduleName string, moduleMap map[string]string, amlData, allData map[string][]float64) ModuleStats {
	// Get genes in this module
	var moduleGenes []string
	for gene, module := range moduleMap {
		if module == moduleName {
			moduleGenes = append(moduleGenes, gene)
		}
	}

	// Get correlation values for both conditions
	amlCorrs := getModuleCorrelations(moduleGenes, amlData)
	allCorrs := getModuleCorrelations(moduleGenes, allData)

	// Calculate t-statistic and p-value manually
	tstat, pval := calculateTTest(amlCorrs, allCorrs)

	return ModuleStats{
		Name: moduleName,
		TStatistic: tstat,
		PValue: pval,
		Size: len(moduleGenes),
	}
}

func getModuleCorrelations(genes []string, expressionData map[string][]float64) []float64 {
	var correlations []float64

	// Get all pairwise correlations
	for i := 0; i < len(genes)-1; i++ {
		for j := i + 1; j < len(genes); j++ {
			gene1 := genes[i]
			gene2 := genes[j]

			// Get expression values
			expr1, ok1 := expressionData[gene1]
			expr2, ok2 := expressionData[gene2]

			if ok1 && ok2 {
				// For ALL data, limit to first 11 samples
				if len(expr1) > 11 {
					expr1 = expr1[:11]
				}
				if len(expr2) > 11 {
					expr2 = expr2[:11]
				}

				// Calculate correlation
				corr := stat.Correlation(expr1, expr2, nil)
				correlations = append(correlations, corr)
			}
		}
	}

	return correlations
}

func calculateTTest(x, y []float64) (tstat, pval float64) {
	meanX := stat.Mean(x, nil)
	meanY := stat.Mean(y, nil)
	varX := stat.Variance(x, nil)
	varY := stat.Variance(y, nil)
	
	nx := float64(len(x))
	ny := float64(len(y))
	
	// Calculate pooled standard error
	se := math.Sqrt((varX/nx) + (varY/ny))
	
	// Calculate t-statistic
	tstat = (meanX - meanY) / se
	
	// Calculate approximate p-value using normal distribution
	z := math.Abs(tstat)
	pval = 2 * (1 - normalCDF(z))
	
	return tstat, pval
}

// normalCDF returns the cumulative distribution function of the standard normal distribution
func normalCDF(x float64) float64 {
	return 0.5 * (1 + math.Erf(x/math.Sqrt(2)))
}

func getUniqueModules(moduleMap map[string]string) map[string]bool {
	modules := make(map[string]bool)
	for _, module := range moduleMap {
		modules[module] = true
	}
	return modules
}
