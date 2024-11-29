package main

import (
	"encoding/csv"
	"math"
	"os"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

// Add the DataWithGenes struct
type DataWithGenes struct {
	Data    *mat.Dense
	GeneIDs []string
}

// Modify ReadData to return gene IDs
func ReadData(filePath string) (*DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1

	// Skip metadata until we find the data table
	for {
		record, err := reader.Read()
		if err != nil {
			return nil, err
		}
		if len(record) > 0 && record[0] == "!dataset_table_begin" {
			break
		}
	}

	// Read header
	headers, err := reader.Read()
	if err != nil {
		return nil, err
	}

	var data [][]string
	var geneIDs []string
	for {
		record, err := reader.Read()
		if err != nil || (len(record) > 0 && record[0] == "!dataset_table_end") {
			break
		}
		geneIDs = append(geneIDs, record[0])
		data = append(data, record)
	}

	numRows := 15923
	numCols := len(headers) - 2
	rawData := mat.NewDense(numRows, numCols, nil)

	for i := 0; i < len(data); i++ {
		for j := 2; j < len(data[i]); j++ {
			value, _ := strconv.ParseFloat(data[i][j], 64)
			rawData.Set(i, j-2, value)
		}
	}

	return &DataWithGenes{Data: rawData, GeneIDs: geneIDs}, nil
}

func ExtractEkerMutants(data *mat.Dense) *mat.Dense {
	// Combine indices 1:12, 25:36, 37:48 as per R code
	indices := make([]int, 36)
	copy(indices[0:12], makeRange(0, 12))
	copy(indices[12:24], makeRange(24, 36))
	copy(indices[24:36], makeRange(36, 48))

	rows, _ := data.Dims()
	ekerData := mat.NewDense(rows, 36, nil)

	for j, idx := range indices {
		col := mat.Col(nil, idx, data)
		for i := 0; i < rows; i++ {
			ekerData.Set(i, j, col[i])
		}
	}

	return ekerData
}

func extractWildTypes(data *mat.Dense) *mat.Dense {
	// Wild types are samples 49:84 in the R code
	rows, _ := data.Dims()
	wildTypeData := mat.NewDense(rows, 36, nil)

	for j := 0; j < 36; j++ {
		col := mat.Col(nil, j+48, data) // 48 is 49-1 for zero-based indexing
		for i := 0; i < rows; i++ {
			wildTypeData.Set(i, j, col[i])
		}
	}

	return wildTypeData
}

// Add ReadGolubData function
func ReadGolubData(filePath string) (*DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1

	// Read header
	headers, err := reader.Read()
	if err != nil {
		return nil, err
	}

	var data [][]string
	var geneIDs []string
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		geneIDs = append(geneIDs, record[0])
		data = append(data, record)
	}

	numRows := len(data)
	numCols := len(headers) - 1
	rawData := mat.NewDense(numRows, numCols, nil)

	for i := 0; i < len(data); i++ {
		for j := 1; j < len(data[i]); j++ {
			value, _ := strconv.ParseFloat(data[i][j], 64)
			rawData.Set(i, j-1, value)
		}
	}

	return &DataWithGenes{Data: rawData, GeneIDs: geneIDs}, nil
}

// Add function to extract ALL samples (similar to ExtractEkerMutants)
func ExtractALLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	allData := mat.NewDense(rows, 27, nil) // ALL has 27 samples

	// Extract ALL samples (columns 0-26)
	for j := 0; j < 27; j++ {
		col := mat.Col(nil, j, data)
		for i := 0; i < rows; i++ {
			allData.Set(i, j, col[i])
		}
	}
	return allData
}

// Add function to extract AML samples (similar to extractWildTypes)
func ExtractAMLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	amlData := mat.NewDense(rows, 11, nil) // AML has 11 samples

	// Extract AML samples (columns 27-37)
	for j := 0; j < 11; j++ {
		col := mat.Col(nil, j+27, data)
		for i := 0; i < rows; i++ {
			amlData.Set(i, j, col[i])
		}
	}
	return amlData
}

// Modify saveToCSV to include gene IDs
func saveToCSV(data *mat.Dense, geneIDs []string, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
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
			return err
		}
	}

	return nil
}

func removeRow(data *mat.Dense, rowIndex int) *mat.Dense {
	rows, cols := data.Dims()
	newData := mat.NewDense(rows-1, cols, nil)

	currentRow := 0
	for i := 0; i < rows; i++ {
		if i != rowIndex {
			row := mat.Row(nil, i, data)
			for j := 0; j < cols; j++ {
				newData.Set(currentRow, j, row[j])
			}
			currentRow++
		}
	}
	return newData
}

func applyLog2(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	logData := mat.NewDense(rows, cols, nil)

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			value := data.At(i, j)
			// Handle zero or negative values
			if value <= 0 {
				logData.Set(i, j, 0)
			} else {
				logData.Set(i, j, math.Log2(value))
			}
		}
	}
	return logData
}

// NormalizeQuantiles applies quantile normalization to the input matrix.
func NormalizeQuantiles(data *mat.Dense) *mat.Dense {
	r, c := data.Dims()
	normData := mat.NewDense(r, c, nil)

	for j := 0; j < c; j++ {
		col := mat.Col(nil, j, data)
		sortFloat(col)
		for i := 0; i < r; i++ {
			normData.Set(i, j, col[i])
		}
	}

	for j := 0; j < c; j++ {
		// Compute the mean for the current column
		col := mat.Col(nil, j, normData)
		mean := meanFloat(col)
		for i := 0; i < r; i++ {
			normData.Set(i, j, normData.At(i, j)-mean)
		}
	}

	return normData
}

// Helper function to sort float slices
func sortFloat(data []float64) {
	for i := 1; i < len(data); i++ {
		key := data[i]
		j := i - 1
		for j >= 0 && data[j] > key {
			data[j+1] = data[j]
			j--
		}
		data[j+1] = key
	}
}

// Helper function to calculate the mean of a float slice
func meanFloat(data []float64) float64 {
	sum := 0.0
	for _, v := range data {
		sum += v
	}
	return sum / float64(len(data))
}

func makeRange(min, max int) []int {
	a := make([]int, max-min)
	for i := range a {
		a[i] = min + i
	}
	return a
}
