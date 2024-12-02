package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestMakeRange(t *testing.T) {
	// Test cases 1 through 4
	for i := 1; i <= 4; i++ {
		// Read input file
		inputPath := filepath.Join("Testing", "makeRange", "In", fmt.Sprintf("input%d.txt", i))
		outputPath := filepath.Join("Testing", "makeRange", "Out", fmt.Sprintf("output%d.txt", i))

		// Read input numbers (min and max)
		min, max, err := readMinMax(inputPath)
		if err != nil {
			t.Errorf("Error reading input file %d: %v", i, err)
			continue
		}

		// Read expected output
		expected, err := readIntSlice(outputPath)
		if err != nil {
			t.Errorf("Error reading output file %d: %v", i, err)
			continue
		}

		// Calculate result
		result := makeRange(min, max)

		// Compare with expected output
		if !sliceEqual(result, expected) {
			t.Errorf("Test case %d failed: got %v, want %v", i, result, expected)
		}
	}
}

func readMinMax(path string) (min, max int, err error) {
	file, err := os.Open(path)
	if err != nil {
		return 0, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if scanner.Scan() {
		min, err = strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return 0, 0, err
		}
	}
	if scanner.Scan() {
		max, err = strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return 0, 0, err
		}
	}
	return min, max, scanner.Err()
}

func readIntSlice(path string) ([]int, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var numbers []int
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		num, err := strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return nil, err
		}
		numbers = append(numbers, num)
	}
	return numbers, scanner.Err()
}

func sliceEqual(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// ReadMatrix reads a matrix from a given file path.
// The first line should contain the row index to remove, followed by the matrix data.
func ReadMatrix(filePath string) (*mat.Dense, int, error) {
	inputData, err := ioutil.ReadFile(filePath)
	if err != nil {
		return nil, 0, err
	}

	lines := strings.Split(string(inputData), "\n")
	if len(lines) < 2 {
		return nil, 0, fmt.Errorf("Input file must contain at least one row of data and one row index")
	}

	// Read the row to remove
	rowToRemove, err := strconv.Atoi(strings.TrimSpace(lines[0]))
	if err != nil {
		return nil, 0, fmt.Errorf("Invalid row index: %v", err)
	}

	// Parse the matrix data
	var dataRows [][]float64
	for _, line := range lines[1:] {
		if line == "" {
			continue
		}
		var row []float64
		for _, val := range strings.Fields(line) {
			num, _ := strconv.ParseFloat(val, 64)
			row = append(row, num)
		}
		dataRows = append(dataRows, row)
	}

	// Create a matrix from the parsed data
	data := mat.NewDense(len(dataRows), len(dataRows[0]), nil)
	for i, row := range dataRows {
		data.SetRow(i, row)
	}

	return data, rowToRemove, nil
}

// MatrixEqual compares two matrices for equality.
func MatrixEqual(a, b *mat.Dense) bool {
	if a == nil || b == nil {
		return a == b
	}
	aRows, aCols := a.Dims()
	bRows, bCols := b.Dims()
	if aRows != bRows || aCols != bCols {
		return false
	}
	for i := 0; i < aRows; i++ {
		for j := 0; j < aCols; j++ {
			if math.Abs(a.At(i, j)-b.At(i, j)) > 1e-9 { // Use a tolerance for floating-point comparison
				return false
			}
		}
	}
	return true
}

func TestRemoveRow(t *testing.T) {
	// Test cases 1 through 4
	for i := 1; i <= 4; i++ {
		// Read input file
		inputPath := filepath.Join("Testing", "removeRow", "In", fmt.Sprintf("input%d.txt", i))
		outputPath := filepath.Join("Testing", "removeRow", "Out", fmt.Sprintf("output%d.txt", i))

		// Read input matrix and row to remove
		matrix, rowToRemove, err := ReadMatrix(inputPath)
		if err != nil {
			t.Errorf("Error reading input file %d: %v", i, err)
			continue
		}

		// Read expected output matrix
		expectedOutput, _, err := ReadMatrix(outputPath)
		if err != nil {
			t.Errorf("Error reading output file %d: %v", i, err)
			continue
		}

		// Remove the specified row
		result := removeRow(matrix, rowToRemove)

		// Compare with expected output
		if !MatrixEqual(result, expectedOutput) {
			t.Errorf("Test case %d failed: matrices are not equal", i)
		}
	}
}
