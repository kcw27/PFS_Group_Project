package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

// DataWithGenes holds both the expression data matrix and gene IDs
type DataWithGenes struct {
	Data    *mat.Dense
	GeneIDs []string
}

// ReadData reads and parses the GDS2901.soft file
func ReadData(filePath string) (*DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1 // Allow variable number of fields
	reader.LazyQuotes = true    // Be more permissive with quotes

	var geneIDs []string
	var dataRows [][]float64
	var dataStarted bool
	var numCols int

	// Read the file line by line
	for {
		record, err := reader.Read()
		if err != nil {
			break // End of file or error
		}

		// Skip empty lines
		if len(record) == 0 {
			continue
		}

		// Look for the start of data
		if record[0] == "!dataset_table_begin" {
			dataStarted = true
			// Read the header in the next iteration
			continue
		}

		// If we found the end of data, break
		if record[0] == "!dataset_table_end" {
			break
		}

		// Skip until we find the start of data
		if !dataStarted {
			continue
		}

		// If this is the first line after data start, it's the header
		if numCols == 0 {
			numCols = len(record) - 2 // Subtract ID and description columns
			continue
		}

		// Process data rows
		if len(record) >= 3 { // Ensure we have at least ID, description, and one data point
			geneIDs = append(geneIDs, record[0])
			rowData := make([]float64, numCols)

			for j := 2; j < len(record) && j-2 < numCols; j++ {
				val, err := strconv.ParseFloat(record[j], 64)
				if err != nil {
					val = 0 // Use 0 for invalid values
				}
				rowData[j-2] = val
			}
			dataRows = append(dataRows, rowData)
		}
	}

	// Verify we have data
	if len(dataRows) == 0 {
		return nil, fmt.Errorf("no valid data found in file")
	}

	// Create matrix from data
	matrix := mat.NewDense(len(dataRows), numCols, nil)
	for i, row := range dataRows {
		matrix.SetRow(i, row)
	}

	return &DataWithGenes{
		Data:    matrix,
		GeneIDs: geneIDs,
	}, nil
}

// ReadGolubData reads and parses the Golub data file
func ReadGolubData(filePath string) (*DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'

	// Read and skip header
	_, err = reader.Read()
	if err != nil {
		return nil, fmt.Errorf("error reading header: %v", err)
	}

	// Read data rows
	var geneIDs []string
	var dataRows [][]float64

	for {
		record, err := reader.Read()
		if err != nil {
			break
		}

		geneIDs = append(geneIDs, record[0])
		rowData := make([]float64, len(record)-1)
		for j := 1; j < len(record); j++ {
			val, err := strconv.ParseFloat(record[j], 64)
			if err != nil {
				return nil, fmt.Errorf("error parsing float: %v", err)
			}
			rowData[j-1] = val
		}
		dataRows = append(dataRows, rowData)
	}

	// Create matrix from data
	numRows, numCols := len(dataRows), len(dataRows[0])
	matrix := mat.NewDense(numRows, numCols, nil)
	for i, row := range dataRows {
		matrix.SetRow(i, row)
	}

	return &DataWithGenes{
		Data:    matrix,
		GeneIDs: geneIDs,
	}, nil
}

// Apply log2 transformation
func applyLog2(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	result := mat.NewDense(rows, cols, nil)
	result.Apply(func(i, j int, v float64) float64 {
		return math.Log2(v + 1)
	}, data)
	return result
}

// NormalizeQuantiles performs quantile normalization
func NormalizeQuantiles(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	result := mat.NewDense(rows, cols, nil)

	// For each column
	for j := 0; j < cols; j++ {
		// Get column data
		col := mat.Col(nil, j, data)

		// Sort values
		sorted := make([]float64, len(col))
		copy(sorted, col)
		sort.Float64s(sorted)

		// Create rank mapping
		rankMap := make(map[float64]float64)
		for i, v := range sorted {
			rankMap[v] = float64(i)
		}

		// Apply normalization
		for i := 0; i < rows; i++ {
			val := data.At(i, j)
			rank := rankMap[val]
			normalizedVal := sorted[int(rank)]
			result.Set(i, j, normalizedVal)
		}
	}

	return result
}

// Extract samples for different conditions
func ExtractEkerSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	ekerCols := makeRange(0, 36)
	result := mat.NewDense(rows, len(ekerCols), nil)
	for j, col := range ekerCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

func ExtractWildSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	wildCols := makeRange(36, 72)
	result := mat.NewDense(rows, len(wildCols), nil)
	for j, col := range wildCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

func ExtractALLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	allCols := makeRange(0, 27)
	result := mat.NewDense(rows, len(allCols), nil)
	for j, col := range allCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

func ExtractAMLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	amlCols := makeRange(27, 38)
	result := mat.NewDense(rows, len(amlCols), nil)
	for j, col := range amlCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

// Helper functions
func removeRow(data *mat.Dense, rowIndex int) *mat.Dense {
	rows, cols := data.Dims()
	if rowIndex < 0 || rowIndex >= rows {
		return data
	}

	newData := mat.NewDense(rows-1, cols, nil)
	for i := 0; i < rowIndex; i++ {
		row := mat.Row(nil, i, data)
		newData.SetRow(i, row)
	}
	for i := rowIndex + 1; i < rows; i++ {
		row := mat.Row(nil, i, data)
		newData.SetRow(i-1, row)
	}
	return newData
}

func makeRange(min, max int) []int {
	a := make([]int, max-min)
	for i := range a {
		a[i] = min + i
	}
	return a
}

func saveToCSV(data *mat.Dense, geneIDs []string, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error creating file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	rows, cols := data.Dims()
	for i := 0; i < rows; i++ {
		row := make([]string, cols+1)
		row[0] = geneIDs[i]
		for j := 0; j < cols; j++ {
			row[j+1] = strconv.FormatFloat(data.At(i, j), 'f', -1, 64)
		}
		if err := writer.Write(row); err != nil {
			return fmt.Errorf("error writing row: %v", err)
		}
	}

	return nil
}
