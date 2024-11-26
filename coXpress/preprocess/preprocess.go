package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

func ReadGolubData(filePath string) (*mat.Dense, error) {
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
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		data = append(data, record)
	}

	numRows := len(data)
	numCols := len(headers) - 1 // Exclude Gene column
	rawData := mat.NewDense(numRows, numCols, nil)

	for i := 0; i < len(data); i++ {
		for j := 1; j < len(data[i]); j++ { // Start from 1 to skip Gene column
			value, _ := strconv.ParseFloat(data[i][j], 64)
			rawData.Set(i, j-1, value)
		}
	}

	return rawData, nil
}

func ExtractALLSamples(data *mat.Dense) *mat.Dense {
	// ALL samples are columns 0-26 (27 samples)
	rows, _ := data.Dims()
	allData := mat.NewDense(rows, 27, nil)

	for j := 0; j < 27; j++ {
		col := mat.Col(nil, j, data)
		for i := 0; i < rows; i++ {
			allData.Set(i, j, col[i])
		}
	}

	return allData
}

func ExtractAMLSamples(data *mat.Dense) *mat.Dense {
	// AML samples are columns 27-37 (11 samples)
	rows, _ := data.Dims()
	amlData := mat.NewDense(rows, 11, nil)

	for j := 0; j < 11; j++ {
		col := mat.Col(nil, j+27, data)
		for i := 0; i < rows; i++ {
			amlData.Set(i, j, col[i])
		}
	}

	return amlData
}

func main() {
	dataFile := "golub.txt"

	rawData, err := ReadGolubData(dataFile)
	if err != nil {
		log.Fatalf("Error reading data: %v", err)
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

	fmt.Println("Data processing complete. Files saved: all_samples.csv, aml_samples.csv")
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
