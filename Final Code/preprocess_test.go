package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

func TestMeanFloat(t *testing.T) {
	// Test cases 1 through 4
	for i := 1; i <= 4; i++ {
		// Read input file
		inputPath := filepath.Join("Testing", "meanFloat", "In", fmt.Sprintf("input%d.txt", i))
		outputPath := filepath.Join("Testing", "meanFloat", "Out", fmt.Sprintf("output%d.txt", i))

		// Read input numbers
		input, err := readFloatSlice(inputPath)
		if err != nil {
			t.Errorf("Error reading input file %d: %v", i, err)
			continue
		}

		// Read expected output
		expected, err := readFloat(outputPath)
		if err != nil {
			t.Errorf("Error reading output file %d: %v", i, err)
			continue
		}

		// Calculate result
		result := meanFloat(input)

		// Compare with expected output
		if result != expected {
			t.Errorf("Test case %d failed: got %v, want %v", i, result, expected)
		}
	}
}

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

func TestSortFloat(t *testing.T) {
	for i := 1; i <= 4; i++ {
		inputFile := "testing/sortFloat/In/input" + strconv.Itoa(i) + ".txt"
		outputFile := "testing/sortFloat/Out/output" + strconv.Itoa(i) + ".txt"

		inputData, err := readFloatSliceFromFile(inputFile)
		if err != nil {
			t.Fatalf("Failed to read input file %s: %v", inputFile, err)
		}

		expectedOutput, err := readFloatSliceFromFile(outputFile)
		if err != nil {
			t.Fatalf("Failed to read output file %s: %v", outputFile, err)
		}

		sortFloat(inputData)

		if !equalFloatSlices(inputData, expectedOutput) {
			t.Errorf("Test %d failed: expected %v, got %v", i, expectedOutput, inputData)
		}
	}
}

// Helper functions for reading test files
func readFloatSlice(path string) ([]float64, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var numbers []float64
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		num, err := strconv.ParseFloat(strings.TrimSpace(scanner.Text()), 64)
		if err != nil {
			return nil, err
		}
		numbers = append(numbers, num)
	}
	return numbers, scanner.Err()
}

func readFloat(path string) (float64, error) {
	file, err := os.Open(path)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if scanner.Scan() {
		return strconv.ParseFloat(strings.TrimSpace(scanner.Text()), 64)
	}
	return 0, scanner.Err()
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

func readFloatSliceFromFile(filename string) ([]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var data []float64
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, ",")
		for _, v := range values {
			if v == "" {
				continue
			}
			num, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
			if err != nil {
				return nil, err
			}
			data = append(data, num)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return data, nil
}

func equalFloatSlices(a, b []float64) bool {
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
