package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

func ReadData(filePath string) (*mat.Dense, error) {
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
	for {
		record, err := reader.Read()
		if err != nil || (len(record) > 0 && record[0] == "!dataset_table_end") {
			break
		}
		data = append(data, record)
	}

	// Create matrix (15923 rows as per R code)
	numRows := 15923
	numCols := len(headers) - 2 // Exclude ID_REF and description columns
	rawData := mat.NewDense(numRows, numCols, nil)

	for i := 0; i < len(data); i++ {
		for j := 2; j < len(data[i]); j++ {
			value, _ := strconv.ParseFloat(data[i][j], 64)
			rawData.Set(i, j-2, value)
		}
	}

	return rawData, nil
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

func makeRange(min, max int) []int {
	a := make([]int, max-min)
	for i := range a {
		a[i] = min + i
	}
	return a
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
}

func saveToCSV(data *mat.Dense, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	rows, cols := data.Dims()
	for i := 0; i < rows; i++ {
		row := make([]string, cols)
		for j := 0; j < cols; j++ {
			row[j] = strconv.FormatFloat(data.At(i, j), 'f', -1, 64)
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	return nil
}
